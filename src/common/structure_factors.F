!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
! **************************************************************************************************
MODULE structure_factors

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE structure_factor_types,          ONLY: structure_factor_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'structure_factors'

   PRIVATE
   PUBLIC :: structure_factor_evaluate, structure_factor_allocate
   PUBLIC :: structure_factor_deallocate, structure_factor_init

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param exp_igr ...
! **************************************************************************************************
   SUBROUTINE structure_factor_init(exp_igr)

      TYPE(structure_factor_type), INTENT(INOUT)         :: exp_igr

      NULLIFY (exp_igr%ex, exp_igr%ey, exp_igr%ez)
      NULLIFY (exp_igr%shell_ex, exp_igr%shell_ey, exp_igr%shell_ez)
      NULLIFY (exp_igr%core_ex, exp_igr%core_ey, exp_igr%core_ez)
      NULLIFY (exp_igr%centre, exp_igr%shell_centre, exp_igr%core_centre)
      NULLIFY (exp_igr%delta, exp_igr%shell_delta, exp_igr%core_delta)

   END SUBROUTINE structure_factor_init

! **************************************************************************************************
!> \brief ...
!> \param exp_igr ...
! **************************************************************************************************
   SUBROUTINE structure_factor_deallocate(exp_igr)

      TYPE(structure_factor_type), INTENT(INOUT)         :: exp_igr

      DEALLOCATE (exp_igr%ex)
      DEALLOCATE (exp_igr%ey)
      DEALLOCATE (exp_igr%ez)
      IF (ASSOCIATED(exp_igr%shell_ex)) THEN
         DEALLOCATE (exp_igr%shell_ex)
         DEALLOCATE (exp_igr%shell_ey)
         DEALLOCATE (exp_igr%shell_ez)
      END IF
      IF (ASSOCIATED(exp_igr%core_ex)) THEN
         DEALLOCATE (exp_igr%core_ex)
         DEALLOCATE (exp_igr%core_ey)
         DEALLOCATE (exp_igr%core_ez)
      END IF
      IF (ASSOCIATED(exp_igr%centre)) THEN
         DEALLOCATE (exp_igr%centre, exp_igr%delta)
      END IF
      IF (ASSOCIATED(exp_igr%shell_centre)) THEN
         DEALLOCATE (exp_igr%shell_centre, exp_igr%shell_delta)
      END IF
      IF (ASSOCIATED(exp_igr%core_centre)) THEN
         DEALLOCATE (exp_igr%core_centre, exp_igr%core_delta)
      END IF

   END SUBROUTINE structure_factor_deallocate

! **************************************************************************************************
!> \brief ...
!> \param bds ...
!> \param nparts ...
!> \param exp_igr ...
!> \param allocate_centre ...
!> \param allocate_shell_e ...
!> \param allocate_shell_centre ...
!> \param nshell ...
! **************************************************************************************************
   SUBROUTINE structure_factor_allocate(bds, nparts, exp_igr, &
                                        allocate_centre, allocate_shell_e, &
                                        allocate_shell_centre, nshell)

      INTEGER, DIMENSION(:, :), INTENT(IN)               :: bds
      INTEGER, INTENT(IN)                                :: nparts
      TYPE(structure_factor_type), INTENT(OUT)           :: exp_igr
      LOGICAL, INTENT(IN), OPTIONAL                      :: allocate_centre, allocate_shell_e, &
                                                            allocate_shell_centre
      INTEGER, INTENT(IN), OPTIONAL                      :: nshell

      ALLOCATE (exp_igr%ex(bds(1, 1):bds(2, 1) + 1, nparts))
      ALLOCATE (exp_igr%ey(bds(1, 2):bds(2, 2) + 1, nparts))
      ALLOCATE (exp_igr%ez(bds(1, 3):bds(2, 3) + 1, nparts))
      NULLIFY (exp_igr%centre, exp_igr%delta)

      exp_igr%lb(1) = LBOUND(exp_igr%ex, 1)
      exp_igr%lb(2) = LBOUND(exp_igr%ey, 1)
      exp_igr%lb(3) = LBOUND(exp_igr%ez, 1)

      IF (PRESENT(allocate_centre)) THEN
         IF (allocate_centre) THEN
            ALLOCATE (exp_igr%centre(3, nparts), exp_igr%delta(3, nparts))
         END IF
      END IF

      IF (PRESENT(allocate_shell_e)) THEN
         IF (allocate_shell_e) THEN
            ALLOCATE (exp_igr%shell_ex(bds(1, 1):bds(2, 1) + 1, nshell))
            ALLOCATE (exp_igr%shell_ey(bds(1, 2):bds(2, 2) + 1, nshell))
            ALLOCATE (exp_igr%shell_ez(bds(1, 3):bds(2, 3) + 1, nshell))
            NULLIFY (exp_igr%shell_centre, exp_igr%shell_delta)

            ALLOCATE (exp_igr%core_ex(bds(1, 1):bds(2, 1) + 1, nshell))
            ALLOCATE (exp_igr%core_ey(bds(1, 2):bds(2, 2) + 1, nshell))
            ALLOCATE (exp_igr%core_ez(bds(1, 3):bds(2, 3) + 1, nshell))
            NULLIFY (exp_igr%core_centre, exp_igr%core_delta)

            IF (PRESENT(allocate_shell_centre)) THEN
               IF (allocate_shell_centre) THEN
                  ALLOCATE (exp_igr%shell_centre(3, nshell), exp_igr%shell_delta(3, nshell))
                  ALLOCATE (exp_igr%core_centre(3, nshell), exp_igr%core_delta(3, nshell))
               END IF
            END IF
         END IF
      ELSE
         NULLIFY (exp_igr%shell_ex, exp_igr%shell_ey, exp_igr%shell_ez)
         NULLIFY (exp_igr%core_ex, exp_igr%core_ey, exp_igr%core_ez)
         NULLIFY (exp_igr%shell_centre, exp_igr%core_centre)
         NULLIFY (exp_igr%shell_delta, exp_igr%core_delta)
      END IF

   END SUBROUTINE structure_factor_allocate

! **************************************************************************************************
!> \brief ...
!> \param delta ...
!> \param lb ...
!> \param ex ...
!> \param ey ...
!> \param ez ...
! **************************************************************************************************
   SUBROUTINE structure_factor_evaluate(delta, lb, ex, ey, ez)

      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: delta
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb
      COMPLEX(KIND=dp), DIMENSION(lb(1):), INTENT(out)   :: ex
      COMPLEX(KIND=dp), DIMENSION(lb(2):), INTENT(out)   :: ey
      COMPLEX(KIND=dp), DIMENSION(lb(3):), INTENT(out)   :: ez

      COMPLEX(KIND=dp)                                   :: fm, fp
      INTEGER                                            :: j, l0, l1, m0, m1, n0, n1
      REAL(KIND=dp)                                      :: vec(3)

      l0 = LBOUND(ex, 1)
      l1 = UBOUND(ex, 1)
      m0 = LBOUND(ey, 1)
      m1 = UBOUND(ey, 1)
      n0 = LBOUND(ez, 1)
      n1 = UBOUND(ez, 1)

      ! delta is in scaled coordinates
      vec(:) = twopi*(delta(:) + 0.5_dp)

      ex(l0) = 1.0_dp
      ey(m0) = 1.0_dp
      ez(n0) = 1.0_dp
      ex(l1) = 1.0_dp
      ey(m1) = 1.0_dp
      ez(n1) = 1.0_dp

      fp = CMPLX(COS(vec(1)), -SIN(vec(1)), KIND=dp)
      fm = CONJG(fp)
      DO j = 1, -l0
         ex(j + l0) = ex(j + l0 - 1)*fp
         ex(-j + l1) = ex(-j + l1 + 1)*fm
      END DO

      fp = CMPLX(COS(vec(2)), -SIN(vec(2)), KIND=dp)
      fm = CONJG(fp)
      DO j = 1, -m0
         ey(j + m0) = ey(j + m0 - 1)*fp
         ey(-j + m1) = ey(-j + m1 + 1)*fm
      END DO

      fp = CMPLX(COS(vec(3)), -SIN(vec(3)), KIND=dp)
      fm = CONJG(fp)
      DO j = 1, -n0
         ez(j + n0) = ez(j + n0 - 1)*fp
         ez(-j + n1) = ez(-j + n1 + 1)*fm
      END DO

   END SUBROUTINE structure_factor_evaluate

END MODULE structure_factors
