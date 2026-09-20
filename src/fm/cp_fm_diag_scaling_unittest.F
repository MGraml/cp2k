!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Eigensolver cost against matrix order and the number of ranks it runs on.
!>        The matrix is the (2,-1) tridiagonal, whose spectrum is known in closed form
!>        (lambda_n = 2 - 2 cos(n pi/(K+1))), so every solve is also checked. Ranks outside
!>        the working subgroup idle, which is what cp_fm_redistribute does, except that here
!>        the matrix is BUILT on the small grid instead of copied onto it.
!>        Rank 0 writes one tab-separated row per (K, P).
!> \author generated for the BSE iterative campaign, 2026-09-19
! **************************************************************************************************
PROGRAM cp_fm_diag_scaling_unittest
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_diag,                      ONLY: FM_DIAG_TYPE_ELPA,&
                                              FM_DIAG_TYPE_SCALAPACK,&
                                              choose_eigv_solver,&
                                              diag_type
   USE cp_fm_elpa,                      ONLY: finalize_elpa_library,&
                                              initialize_elpa_library
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE kinds,                           ONLY: dp,&
                                              int_4,&
                                              int_8
   USE machine,                         ONLY: default_output_unit,&
                                              m_memory,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init

!$ USE OMP_LIB, ONLY: omp_get_max_threads
#include "../base/base_uses.f90"

   IMPLICIT NONE

   INTEGER, PARAMETER, DIMENSION(6) :: k_default = [256, 512, 1024, 2048, 4096, 8192]
   INTEGER, PARAMETER :: max_rows_per_rank = 2048   ! skip combinations slower than this
   REAL(KIND=dp), PARAMETER :: tol = 1.0E-8_dp

   INTEGER :: ik, n_k, k, p, npe, mepos, io_unit, ngroups, i, j, ir, jc, nrl, ncl, nprow, npcol, &
              nrb, ncb, nthr, grid_r, grid_c, n_over
   INTEGER, ALLOCATABLE, DIMENSION(:) :: k_list
   INTEGER(KIND=int_8) :: loc_elem, glob_elem, panel_col, panel_row
   CHARACTER(LEN=16) :: arg
   CHARACTER(LEN=256) :: karg
   LOGICAL :: use_elpa
   INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices, group_distribution, group_partition
   REAL(KIND=dp) :: t0, dt, lam1_ref, err
   INTEGER(KIND=int_8) :: mem0, mem1
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigval
   TYPE(cp_blacs_env_type), POINTER :: blacs_env
   TYPE(cp_fm_struct_type), POINTER :: fm_struct
   TYPE(cp_fm_type) :: fm_a, fm_vec
   TYPE(mp_para_env_type), POINTER :: para_env, sub_env
   TYPE(cp_logger_type), POINTER :: logger

   NULLIFY (para_env, sub_env, blacs_env, fm_struct, group_distribution, group_partition, logger)
   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   npe = para_env%num_pe
   mepos = para_env%mepos
   nthr = 1
!$ nthr = omp_get_max_threads()
   io_unit = default_output_unit
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io_unit, &
                         close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)

   ! backend from the command line: default ScaLAPACK, "elpa" for ELPA
   arg = ""
   IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, arg)
   use_elpa = (TRIM(arg) == "elpa")
   IF (use_elpa) THEN
      diag_type = FM_DIAG_TYPE_ELPA
      CALL initialize_elpa_library()
   ELSE
      diag_type = FM_DIAG_TYPE_SCALAPACK
   END IF

   ! optional forced process grid: "<backend> <nprow> <npcol>" pins the BLACS layout and runs the
   ! single rank count nprow*npcol, which is how a skewed inherited grid is reproduced cheaply
   grid_r = 0
   grid_c = 0
   IF (COMMAND_ARGUMENT_COUNT() >= 3) THEN
      CALL GET_COMMAND_ARGUMENT(2, arg)
      READ (arg, *) grid_r
      CALL GET_COMMAND_ARGUMENT(3, arg)
      READ (arg, *) grid_c
      IF (grid_r*grid_c > npe) THEN
         CALL cp_abort(__LOCATION__, &
                       "the forced process grid asks for more ranks than the job has")
      END IF
   END IF

   ! optional 4th argument: a comma-separated list of matrix orders, replacing the built-in
   ! ladder, so one production-sized order can be asked for without editing the source
   karg = ""
   IF (COMMAND_ARGUMENT_COUNT() >= 4) CALL GET_COMMAND_ARGUMENT(4, karg)
   IF (LEN_TRIM(karg) > 0) THEN
      n_k = 1
      DO i = 1, LEN_TRIM(karg)
         IF (karg(i:i) == ",") THEN
            n_k = n_k + 1
            karg(i:i) = " "
         END IF
      END DO
      ALLOCATE (k_list(n_k))
      READ (karg, *) k_list
   ELSE
      n_k = SIZE(k_default)
      ALLOCATE (k_list(n_k))
      k_list = k_default
   END IF

   IF (mepos == 0) THEN
      WRITE (io_unit, '(A)') "FMSCAL"//CHAR(9)//"K"//CHAR(9)//"P"//CHAR(9)//"npe"//CHAR(9)// &
         "nthr"//CHAR(9)//"nprow"//CHAR(9)//"npcol"//CHAR(9)//"time_s"//CHAR(9)//"mem_mb"// &
         CHAR(9)//"err"
      WRITE (io_unit, '(A)') "FMGEOM"//CHAR(9)//"K"//CHAR(9)//"P"//CHAR(9)//"nprow"//CHAR(9)// &
         "npcol"//CHAR(9)//"nrb"//CHAR(9)//"ncb"//CHAR(9)//"nrl"//CHAR(9)//"ncl"//CHAR(9)// &
         "loc_elem"//CHAR(9)//"panel_col"//CHAR(9)//"panel_row"//CHAR(9)//"rss_mb"//CHAR(9)// &
         "n_over_i4"
   END IF

   DO ik = 1, n_k
      k = k_list(ik)
      p = 1
      IF (grid_r > 0) p = grid_r*grid_c
      DO WHILE (p <= npe)
         IF (k/p <= max_rows_per_rank .OR. LEN_TRIM(karg) > 0) THEN
            ! first p ranks form the working group; the rest idle, as in cp_fm_redistribute
            IF (p < npe) THEN
               ALLOCATE (group_distribution(0:npe - 1), group_partition(0:1))
               group_partition = [p, npe - p]
               ALLOCATE (sub_env)
               CALL sub_env%from_split(comm=para_env, ngroups=ngroups, &
                                       group_distribution=group_distribution, &
                                       n_subgroups=2, group_partition=group_partition)
            END IF

            IF (p == npe .OR. group_distribution(mepos) == 0) THEN
               IF (p == npe) THEN
                  IF (grid_r > 0) THEN
                     CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env, &
                                              grid_2d=[grid_r, grid_c])
                  ELSE
                     CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
                  END IF
                  CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                           nrow_global=k, ncol_global=k)
               ELSE
                  IF (grid_r > 0) THEN
                     CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=sub_env, &
                                              grid_2d=[grid_r, grid_c])
                  ELSE
                     CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=sub_env)
                  END IF
                  CALL cp_fm_struct_create(fm_struct, para_env=sub_env, context=blacs_env, &
                                           nrow_global=k, ncol_global=k)
               END IF
               CALL cp_fm_create(fm_a, fm_struct, name="tridiag")
               CALL cp_fm_create(fm_vec, fm_struct, name="vectors")
               CALL cp_fm_struct_release(fm_struct)

               ! local fill of the (2,-1) tridiagonal
               CALL cp_fm_get_info(fm_a, nrow_local=nrl, ncol_local=ncl, &
                                   row_indices=row_indices, col_indices=col_indices, &
                                   nrow_block=nrb, ncol_block=ncb)
               fm_a%local_data(:, :) = 0.0_dp
               DO jc = 1, ncl
                  j = col_indices(jc)
                  DO ir = 1, nrl
                     i = row_indices(ir)
                     IF (i == j) fm_a%local_data(ir, jc) = 2.0_dp
                     IF (ABS(i - j) == 1) fm_a%local_data(ir, jc) = -1.0_dp
                  END DO
               END DO

               ALLOCATE (eigval(k))
               eigval(:) = 0.0_dp

               ! the sizes an INTEGER(4) has to hold on this geometry, worst rank, before the
               ! solve is entered: at n_ov = 275176 on a 2 x 240 grid a row panel is 3.8E+10
               ! elements, so a report that only appears after the solve reports nothing
               IF (p == npe) THEN
                  CALL para_env%max(nrl)
                  CALL para_env%max(ncl)
               ELSE
                  CALL sub_env%max(nrl)
                  CALL sub_env%max(ncl)
               END IF
               loc_elem = INT(nrl, int_8)*INT(ncl, int_8)
               glob_elem = INT(k, int_8)*INT(k, int_8)
               panel_col = INT(k, int_8)*INT(ncl, int_8)
               panel_row = INT(nrl, int_8)*INT(k, int_8)
               n_over = 0
               IF (loc_elem > INT(HUGE(0_int_4), int_8)) n_over = n_over + 1
               IF (glob_elem > INT(HUGE(0_int_4), int_8)) n_over = n_over + 1
               IF (panel_col > INT(HUGE(0_int_4), int_8)) n_over = n_over + 1
               IF (panel_row > INT(HUGE(0_int_4), int_8)) n_over = n_over + 1
               IF (8_int_8*loc_elem > INT(HUGE(0_int_4), int_8)) n_over = n_over + 1
               nprow = SIZE(blacs_env%blacs2mpi, 1)
               npcol = SIZE(blacs_env%blacs2mpi, 2)
               CALL m_memory(mem0)
               IF (mepos == 0) THEN
                  WRITE (io_unit, '(A,7(A1,I0),4(A1,I0),A1,F10.1,A1,I0)') "FMGEOM", &
                     CHAR(9), k, CHAR(9), p, CHAR(9), nprow, CHAR(9), npcol, CHAR(9), nrb, &
                     CHAR(9), ncb, CHAR(9), nrl, CHAR(9), ncl, &
                     CHAR(9), loc_elem, CHAR(9), panel_col, CHAR(9), panel_row, &
                     CHAR(9), REAL(mem0, dp)/1048576.0_dp, CHAR(9), n_over
                  FLUSH (io_unit)
               END IF

               IF (p == npe) THEN
                  CALL para_env%sync()
               ELSE
                  CALL sub_env%sync()
               END IF
               t0 = m_walltime()
               CALL choose_eigv_solver(fm_a, fm_vec, eigval)
               dt = m_walltime() - t0
               CALL m_memory(mem1)

               lam1_ref = 2.0_dp - 2.0_dp*COS(ACOS(-1.0_dp)/REAL(k + 1, dp))
               err = ABS(eigval(1) - lam1_ref)
               IF (mepos == 0) THEN
                  WRITE (io_unit, '(A,6(A1,I0),A1,ES12.5,A1,F10.1,A1,ES12.5)') "FMSCAL", &
                     CHAR(9), k, CHAR(9), p, CHAR(9), npe, CHAR(9), nthr, CHAR(9), nprow, &
                     CHAR(9), npcol, &
                     CHAR(9), dt, CHAR(9), REAL(mem1 - mem0, dp)/1048576.0_dp, CHAR(9), err
               END IF
               DEALLOCATE (eigval)
               CALL cp_fm_release(fm_a)
               CALL cp_fm_release(fm_vec)
               CALL cp_blacs_env_release(blacs_env)
               IF (err > tol*REAL(k, dp)) THEN
                  CALL cp_abort(__LOCATION__, &
                                "lowest eigenvalue of the tridiagonal test matrix is wrong")
               END IF
            END IF

            IF (p < npe) THEN
               CALL sub_env%free()
               DEALLOCATE (sub_env)
               NULLIFY (sub_env)
               DEALLOCATE (group_distribution, group_partition)
               NULLIFY (group_distribution, group_partition)
            END IF
            CALL para_env%sync()
         END IF
         IF (grid_r > 0) THEN
            EXIT      ! a forced grid fixes the rank count, so there is nothing to sweep
         ELSE IF (p < npe .AND. p*2 > npe) THEN
            p = npe   ! npe is not a power of two: sample it anyway, it is the default geometry
         ELSE
            p = p*2
         END IF
      END DO
   END DO

   DEALLOCATE (k_list)
   IF (use_elpa) CALL finalize_elpa_library()
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   IF (mepos == 0) WRITE (io_unit, '(A)') "FMSCAL done"
   CALL mp_world_finalize()

END PROGRAM cp_fm_diag_scaling_unittest
