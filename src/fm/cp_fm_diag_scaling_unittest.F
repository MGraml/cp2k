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
   USE cp_blacs_env, ONLY: cp_blacs_env_create, cp_blacs_env_release, cp_blacs_env_type
   USE cp_fm_diag, ONLY: FM_DIAG_TYPE_ELPA, FM_DIAG_TYPE_SCALAPACK, choose_eigv_solver, diag_type
   USE cp_fm_elpa, ONLY: finalize_elpa_library, initialize_elpa_library
   USE cp_fm_struct, ONLY: cp_fm_struct_create, cp_fm_struct_release, cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, cp_fm_get_info, cp_fm_release, cp_fm_type
   USE cp_log_handling, ONLY: cp_add_default_logger, cp_logger_create, cp_logger_release, &
                              cp_logger_type, cp_rm_default_logger
   USE kinds, ONLY: dp, int_8
   USE machine, ONLY: default_output_unit, m_memory, m_walltime
   USE message_passing, ONLY: mp_para_env_type, mp_world_finalize, mp_world_init
#include "../base/base_uses.f90"

   IMPLICIT NONE

   INTEGER, PARAMETER :: n_k = 6
   INTEGER, PARAMETER, DIMENSION(n_k) :: k_list = [256, 512, 1024, 2048, 4096, 8192]
   INTEGER, PARAMETER :: max_rows_per_rank = 2048   ! skip combinations slower than this
   REAL(KIND=dp), PARAMETER :: tol = 1.0E-8_dp

   INTEGER :: ik, ip, k, p, npe, mepos, io_unit, ngroups, i, j, ir, jc, nrl, ncl, nprow, npcol
   CHARACTER(LEN=16) :: arg
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

   IF (mepos == 0) THEN
      WRITE (io_unit, '(A)') "FMSCAL"//CHAR(9)//"K"//CHAR(9)//"P"//CHAR(9)//"npe"//CHAR(9)// &
         "nprow"//CHAR(9)//"npcol"//CHAR(9)//"time_s"//CHAR(9)//"mem_mb"//CHAR(9)//"err"
   END IF

   DO ik = 1, n_k
      k = k_list(ik)
      p = 1
      DO WHILE (p <= npe)
         IF (k/p <= max_rows_per_rank) THEN
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
                  CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
                  CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                           nrow_global=k, ncol_global=k)
               ELSE
                  CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=sub_env)
                  CALL cp_fm_struct_create(fm_struct, para_env=sub_env, context=blacs_env, &
                                           nrow_global=k, ncol_global=k)
               END IF
               CALL cp_fm_create(fm_a, fm_struct, name="tridiag")
               CALL cp_fm_create(fm_vec, fm_struct, name="vectors")
               CALL cp_fm_struct_release(fm_struct)

               ! local fill of the (2,-1) tridiagonal
               CALL cp_fm_get_info(fm_a, nrow_local=nrl, ncol_local=ncl, &
                                   row_indices=row_indices, col_indices=col_indices, &
                                   nrow_block=nprow, ncol_block=npcol)
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
               CALL m_memory(mem0)
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
               nprow = SIZE(blacs_env%blacs2mpi, 1)
               npcol = SIZE(blacs_env%blacs2mpi, 2)
               IF (mepos == 0) THEN
                  WRITE (io_unit, '(A,5(A1,I0),A1,ES12.5,A1,F10.1,A1,ES12.5)') "FMSCAL", &
                     CHAR(9), k, CHAR(9), p, CHAR(9), npe, CHAR(9), nprow, CHAR(9), npcol, &
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
         p = p*2
      END DO
   END DO

   IF (use_elpa) CALL finalize_elpa_library()
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   IF (mepos == 0) WRITE (io_unit, '(A)') "FMSCAL done"
   CALL mp_world_finalize()

END PROGRAM cp_fm_diag_scaling_unittest
