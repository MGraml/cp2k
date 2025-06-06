!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines using linear scaling chebyshev methods
!> \par History
!>       2012.10 created [Jinwoong Cha]
!> \author Jinwoong Cha
! **************************************************************************************************
MODULE dm_ls_chebyshev
   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_get_info, dbcsr_get_occupation, dbcsr_multiply, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_trace
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE dm_ls_scf_qs,                    ONLY: write_matrix_to_cube
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: pi
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_chebyshev'

   PUBLIC :: compute_chebyshev

CONTAINS

! **************************************************************************************************
!> \brief compute chebyshev polynomials up to order n for a given value of x
!> \param value ...
!> \param x ...
!> \param n ...
!> \par History
!>       2012.11 created [Jinwoong Cha]
!> \author Jinwoong Cha
! **************************************************************************************************
   SUBROUTINE chebyshev_poly(value, x, n)
      REAL(KIND=dp), INTENT(OUT)                         :: value
      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: n

!polynomial values
!number of chev polynomials

      value = COS((n - 1)*ACOS(x))

   END SUBROUTINE chebyshev_poly

! **************************************************************************************************
!> \brief kernel for chebyshev polynomials expansion (Jackson kernel)
!> \param value ...
!> \param n ...
!> \param nc ...
!> \par History
!>       2012.11 created [Jinwoong Cha]
!> \author Jinwoong Cha
! **************************************************************************************************
   SUBROUTINE kernel(value, n, nc)
      REAL(KIND=dp), INTENT(OUT)                         :: value
      INTEGER, INTENT(IN)                                :: n, nc

!kernel at n
!n-1 order of chebyshev polynomials
!number of total chebyshev polynomials
!Kernel define

      value = 1.0_dp/(nc + 1.0_dp)*((nc - (n - 1) + 1.0_dp)* &
                                    COS(pi*(n - 1)/(nc + 1.0_dp)) + SIN(pi*(n - 1)/(nc + 1.0_dp))*1.0_dp/TAN(pi/(nc + 1.0_dp)))

   END SUBROUTINE kernel

! **************************************************************************************************
!> \brief compute properties based on chebyshev expansion
!> \param qs_env ...
!> \param ls_scf_env ...
!> \par History
!>       2012.10 created [Jinwoong Cha]
!> \author Jinwoong Cha
! **************************************************************************************************
   SUBROUTINE compute_chebyshev(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_chebyshev'
      REAL(KIND=dp), PARAMETER                           :: scale_evals = 1.01_dp

      CHARACTER(LEN=30)                                  :: middle_name
      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: handle, icheb, igrid, iinte, ispin, &
                                                            iwindow, n_gridpoint_dos, ncheb, &
                                                            ninte, Nrows, nwindow, unit_cube, &
                                                            unit_dos, unit_nr
      LOGICAL                                            :: converged, write_cubes
      REAL(KIND=dp) :: chev_T, chev_T_dos, dummy1, final, frob_matrix, initial, interval_a, &
         interval_b, max_ev, min_ev, occ, orbital_occ, summa, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chev_E, chev_Es_dos, dos, dummy2, ev1, &
                                                            ev2, kernel_g, mu, sev1, sev2, trace_dm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aitchev_T, E_inte, gdensity, sqrt_vec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_r
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_dummy1, matrix_F, matrix_tmp1, &
                                                            matrix_tmp2, matrix_tmp3
      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: matrix_dummy2

      IF (.NOT. ls_scf_env%chebyshev%compute_chebyshev) RETURN

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      ncheb = ls_scf_env%chebyshev%n_chebyshev
      ninte = 2*ncheb
      n_gridpoint_dos = ls_scf_env%chebyshev%n_gridpoint_dos

      write_cubes = BTEST(cp_print_key_should_output(logger%iter_info, ls_scf_env%chebyshev%print_key_cube), cp_p_file)
      IF (write_cubes) THEN
         IF (ASSOCIATED(ls_scf_env%chebyshev%min_energy)) DEALLOCATE (ls_scf_env%chebyshev%min_energy)
         CALL section_vals_val_get(ls_scf_env%chebyshev%print_key_cube, "MIN_ENERGY", r_vals=tmp_r)
         ALLOCATE (ls_scf_env%chebyshev%min_energy(SIZE(tmp_r)))
         ls_scf_env%chebyshev%min_energy = tmp_r

         IF (ASSOCIATED(ls_scf_env%chebyshev%max_energy)) DEALLOCATE (ls_scf_env%chebyshev%max_energy)
         CALL section_vals_val_get(ls_scf_env%chebyshev%print_key_cube, "MAX_ENERGY", r_vals=tmp_r)
         ALLOCATE (ls_scf_env%chebyshev%max_energy(SIZE(tmp_r)))
         ls_scf_env%chebyshev%max_energy = tmp_r

         nwindow = SIZE(ls_scf_env%chebyshev%min_energy)
      ELSE
         nwindow = 0
      END IF

      ALLOCATE (ev1(1:nwindow))
      ALLOCATE (ev2(1:nwindow))
      ALLOCATE (sev1(1:nwindow))
      ALLOCATE (sev2(1:nwindow))
      ALLOCATE (trace_dm(1:nwindow))
      ALLOCATE (matrix_dummy2(1:nwindow))

      DO iwindow = 1, nwindow
         ev1(iwindow) = ls_scf_env%chebyshev%min_energy(iwindow)
         ev2(iwindow) = ls_scf_env%chebyshev%max_energy(iwindow)
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A)') "STARTING CHEBYSHEV CALCULATION"
      END IF

      ! create 3 temporary matrices
      CALL dbcsr_create(matrix_tmp1, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp2, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp3, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_F, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_dummy1, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)

      DO iwindow = 1, nwindow
         CALL dbcsr_create(matrix_dummy2(iwindow), template=ls_scf_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      DO ispin = 1, SIZE(ls_scf_env%matrix_ks)
         ! create matrix_F=inv(sqrt(S))*H*inv(sqrt(S))
         CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                             0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                             0.0_dp, matrix_F, filter_eps=ls_scf_env%eps_filter)

         ! find largest and smallest eigenvalues
         CALL arnoldi_extremal(matrix_F, max_ev, min_ev, converged=converged, max_iter=ls_scf_env%max_iter_lanczos, &
                               threshold=ls_scf_env%eps_lanczos) !Lanczos algorithm to calculate eigenvalue
         IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,2F16.8,A,L2)') &
            "smallest largest eigenvalue", min_ev, max_ev, " converged ", converged
         IF (nwindow > 0) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,1000F16.8)') "requested interval-min_energy", ev1(:)
            IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,1000F16.8)') "requested interval-max_energy", ev2(:)
         END IF
         interval_a = (max_ev - min_ev)*scale_evals/2
         interval_b = (max_ev + min_ev)/2

         sev1(:) = (ev1(:) - interval_b)/interval_a !scaled ev1 vector
         sev2(:) = (ev2(:) - interval_b)/interval_a !scaled ev2 vector

         !chebyshev domain,pi*sqrt(1-x^2) vector construction and chebyshev polynomials for integration (for g(E))
         ALLOCATE (E_inte(1:ninte + 1, 1:nwindow))
         ALLOCATE (sqrt_vec(1:ninte + 1, 1:nwindow))

         DO iwindow = 1, nwindow
            DO iinte = 1, ninte + 1
               E_inte(iinte, iwindow) = sev1(iwindow) + ((sev2(iwindow) - sev1(iwindow))/ninte)*(iinte - 1)
               sqrt_vec(iinte, iwindow) = pi*SQRT(1.0_dp - E_inte(iinte, iwindow)*E_inte(iinte, iwindow))
            END DO
         END DO

         !integral.. (identical to the coefficient for g(E))

         ALLOCATE (aitchev_T(1:ncheb, 1:nwindow)) !after intergral. =>ainte

         DO iwindow = 1, nwindow
            DO icheb = 1, ncheb
               CALL chebyshev_poly(initial, E_inte(1, iwindow), icheb)
               CALL chebyshev_poly(final, E_inte(1, iwindow), icheb)
          summa = (sev2(iwindow) - sev1(iwindow))/(2.0_dp*ninte)*(initial/sqrt_vec(1, iwindow) + final/sqrt_vec(ninte + 1, iwindow))
               DO iinte = 2, ninte
                  CALL chebyshev_poly(chev_T, E_inte(iinte, iwindow), icheb)
                  summa = summa + ((sev2(iwindow) - sev1(iwindow))/ninte)*(chev_T/sqrt_vec(iinte, iwindow))
               END DO
               aitchev_T(icheb, iwindow) = summa
               summa = 0
            END DO
         END DO

         ! scale the matrix to get evals in the interval -1,1
         CALL dbcsr_add_on_diag(matrix_F, -interval_b)
         CALL dbcsr_scale(matrix_F, 1/interval_a)

         ! compute chebyshev matrix recursion
         CALL dbcsr_get_info(matrix=matrix_F, nfullrows_total=Nrows) !get information about a matrix
         CALL dbcsr_set(matrix_dummy1, 0.0_dp) !empty matrix creation(for density matrix)

         DO iwindow = 1, nwindow
            CALL dbcsr_set(matrix_dummy2(iwindow), 0.0_dp) !empty matrix creation(for density matrix)
         END DO

         ALLOCATE (mu(1:ncheb))
         ALLOCATE (kernel_g(1:ncheb))
         CALL kernel(kernel_g(1), 1, ncheb)
         CALL kernel(kernel_g(2), 2, ncheb)

         CALL dbcsr_set(matrix_tmp1, 0.0_dp) !matrix creation
         CALL dbcsr_add_on_diag(matrix_tmp1, 1.0_dp) !add a only number to diagonal elements
         CALL dbcsr_trace(matrix_tmp1, trace=mu(1))
         CALL dbcsr_copy(matrix_tmp2, matrix_F) !make matrix_tmp2 = matrix_F
         CALL dbcsr_trace(matrix_tmp2, trace=mu(2))

         DO iwindow = 1, nwindow
            CALL dbcsr_copy(matrix_dummy1, matrix_tmp1)
            CALL dbcsr_copy(matrix_dummy2(iwindow), matrix_tmp2) !matrix_dummy2=
            CALL dbcsr_scale(matrix_dummy1, kernel_g(1)*aitchev_T(1, iwindow)) !first term of chebyshev poly(matrix)
            CALL dbcsr_scale(matrix_dummy2(iwindow), 2.0_dp*kernel_g(2)*aitchev_T(2, iwindow)) !second term of chebyshev poly(matrix)

            CALL dbcsr_add(matrix_dummy2(iwindow), matrix_dummy1, 1.0_dp, 1.0_dp)
         END DO

         DO icheb = 2, ncheb - 1
            t1 = m_walltime()
            CALL dbcsr_multiply("N", "N", 2.0_dp, matrix_F, matrix_tmp2, &
                                -1.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter) !matrix multiplication(Recursion)
            CALL dbcsr_copy(matrix_tmp3, matrix_tmp1)
            CALL dbcsr_copy(matrix_tmp1, matrix_tmp2)
            CALL dbcsr_copy(matrix_tmp2, matrix_tmp3)
            CALL dbcsr_trace(matrix_tmp2, trace=mu(icheb + 1)) !icheb+1 th coefficient
            CALL kernel(kernel_g(icheb + 1), icheb + 1, ncheb)

            DO iwindow = 1, nwindow

               CALL dbcsr_copy(matrix_dummy1, matrix_tmp2)
               CALL dbcsr_scale(matrix_dummy1, 2.0_dp*kernel_g(icheb + 1)*aitchev_T(icheb + 1, iwindow)) !second term of chebyshev poly(matrix)
               CALL dbcsr_add(matrix_dummy2(iwindow), matrix_dummy1, 1.0_dp, 1.0_dp)
               CALL dbcsr_trace(matrix_dummy2(iwindow), trace=trace_dm(iwindow)) !icheb+1 th coefficient

            END DO

            occ = dbcsr_get_occupation(matrix_tmp1)
            t2 = m_walltime()
            IF (unit_nr > 0 .AND. MOD(icheb, 20) == 0) THEN
               CALL m_flush(unit_nr)
               IF (nwindow > 0) THEN
                  WRITE (unit_nr, '(T2,A,I5,1X,A,1X,F8.3,1X,A,1X,F8.6,1X,A,1X,1000F16.8)') &
                     "Iter.", icheb, "time=", t2 - t1, "occ=", occ, "traces=", trace_dm(:)
               ELSE
                  WRITE (unit_nr, '(T2,A,I5,1X,A,1X,F8.3,1X,A,1X,F8.6)') &
                     "Iter.", icheb, "time=", t2 - t1, "occ=", occ
               END IF
            END IF
         END DO

         DO iwindow = 1, nwindow
            IF (SIZE(ls_scf_env%matrix_ks) == 1) THEN
               orbital_occ = 2.0_dp
            ELSE
               orbital_occ = 1.0_dp
            END IF
            CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, matrix_dummy2(iwindow), &
                                0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "N", orbital_occ, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter)
            CALL dbcsr_copy(matrix_dummy2(iwindow), matrix_tmp2)

            ! look at the difference with the density matrix from the ls routines
            IF (.FALSE.) THEN
               CALL dbcsr_copy(matrix_tmp1, matrix_tmp2)
               CALL dbcsr_add(matrix_tmp1, ls_scf_env%matrix_p(ispin), 1.0_dp, -1.0_dp) !comparison
               frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
               IF (unit_nr > 0) WRITE (unit_nr, *) "Difference between Chebyshev DM and LS DM", frob_matrix
            END IF
         END DO

         write_cubes = BTEST(cp_print_key_should_output(logger%iter_info, &
                                                        ls_scf_env%chebyshev%print_key_cube), cp_p_file)
         IF (write_cubes) THEN
            DO iwindow = 1, nwindow
               WRITE (middle_name, "(A,I0)") "E_DENSITY_WINDOW_", iwindow
               WRITE (title, "(A,1X,F16.8,1X,A,1X,F16.8)") "Energy range : ", ev1(iwindow), "to", ev2(iwindow)
               unit_cube = cp_print_key_unit_nr(logger, ls_scf_env%chebyshev%print_key_cube, &
                                                "", extension=".cube", & !added 01/22/2012
                                                middle_name=TRIM(middle_name), log_filename=.FALSE.)
               CALL write_matrix_to_cube(qs_env, ls_scf_env, matrix_dummy2(iwindow), unit_cube, title, &
                                         section_get_ivals(ls_scf_env%chebyshev%print_key_cube, "STRIDE"))
               CALL cp_print_key_finished_output(unit_cube, logger, ls_scf_env%chebyshev%print_key_cube, "")
            END DO
         END IF

      END DO

      ! Chebyshev expansion with calculated coefficient
      ! grid construction and rescaling (by J)
      unit_dos = cp_print_key_unit_nr(logger, ls_scf_env%chebyshev%print_key_dos, "", extension=".xy", &
                                      middle_name="DOS", log_filename=.FALSE.)

      IF (unit_dos > 0) THEN
         ALLOCATE (dos(1:n_gridpoint_dos))
         ALLOCATE (gdensity(1:n_gridpoint_dos, 1:nwindow))
         ALLOCATE (chev_E(1:n_gridpoint_dos))
         ALLOCATE (chev_Es_dos(1:n_gridpoint_dos))
         ALLOCATE (dummy2(1:nwindow))
         DO igrid = 1, n_gridpoint_dos
            chev_E(igrid) = min_ev + (igrid - 1)*(max_ev - min_ev)/(n_gridpoint_dos - 1)
            chev_Es_dos(igrid) = (chev_E(igrid) - interval_b)/interval_a
         END DO
         DO igrid = 1, n_gridpoint_dos
            dummy1 = 0.0_dp !summation of polynomials
            dummy2(:) = 0.0_dp !summation of polynomials
            DO icheb = 2, ncheb
               CALL chebyshev_poly(chev_T_dos, chev_Es_dos(igrid), icheb)
               dummy1 = dummy1 + kernel_g(icheb)*mu(icheb)*chev_T_dos
               DO iwindow = 1, nwindow
                  dummy2(iwindow) = dummy2(iwindow) + kernel_g(icheb)*aitchev_T(icheb, iwindow)*chev_T_dos
               END DO
            END DO
            dos(igrid) = 1.0_dp/(interval_a*Nrows* &
                                 (pi*SQRT(1.0_dp - chev_Es_dos(igrid)*chev_Es_dos(igrid))))*(kernel_g(1)*mu(1) + 2.0_dp*dummy1)
            DO iwindow = 1, nwindow
               gdensity(igrid, iwindow) = kernel_g(1)*aitchev_T(1, iwindow) + 2.0_dp*dummy2(iwindow)
            END DO
            WRITE (unit_dos, '(1000F16.8)') chev_E(igrid), dos(igrid), gdensity(igrid, :)
         END DO
         DEALLOCATE (chev_Es_dos, chev_E, dos, gdensity)
      END IF
      CALL cp_print_key_finished_output(unit_dos, logger, ls_scf_env%chebyshev%print_key_dos, "")

      ! free the matrices
      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)
      CALL dbcsr_release(matrix_tmp3)
      CALL dbcsr_release(matrix_F)
      CALL dbcsr_release(matrix_dummy1)

      DO iwindow = 1, nwindow
         CALL dbcsr_release(matrix_dummy2(iwindow))
      END DO

      DEALLOCATE (ev1, ev2, sev1, sev2, matrix_dummy2)

      !Need deallocation
      DEALLOCATE (mu, kernel_g, aitchev_T, E_inte, sqrt_vec)

      IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') "ENDING CHEBYSHEV CALCULATION"

      CALL timestop(handle)

   END SUBROUTINE compute_chebyshev

END MODULE dm_ls_chebyshev
