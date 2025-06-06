!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for the functional calculations
!> \par History
!>      JGH (20.02.2001) : Added setup routine
!>      JGH (26.02.2003) : OpenMP enabled
!> \author JGH (15.02.2002)
! **************************************************************************************************
MODULE xc_functionals_utilities

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   REAL(KIND=dp), PARAMETER :: rsfac = 0.6203504908994000166680065_dp ! (4*pi/3)^(-1/3)
   REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                               f23 = 2.0_dp*f13, &
                               f43 = 4.0_dp*f13, &
                               f53 = 5.0_dp*f13
   REAL(KIND=dp), PARAMETER :: fxfac = 1.923661050931536319759455_dp ! 1/(2^(4/3) - 2)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_functionals_utilities'

   PUBLIC :: set_util, calc_rs, calc_fx, &
             calc_wave_vector, calc_z, calc_rs_pw, calc_srs_pw

   INTERFACE calc_fx
      MODULE PROCEDURE calc_fx_array, calc_fx_single
   END INTERFACE

   INTERFACE calc_rs
      MODULE PROCEDURE calc_rs_array, calc_rs_single
   END INTERFACE

   REAL(KIND=dp), SAVE :: eps_rho

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE set_util(cutoff)

      REAL(KIND=dp)                                      :: cutoff

      eps_rho = cutoff

   END SUBROUTINE set_util

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param rs ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE calc_rs_single(rho, rs)

!   rs parameter : f*rho**(-1/3)

      REAL(KIND=dp), INTENT(IN)                          :: rho
      REAL(KIND=dp), INTENT(OUT)                         :: rs

      IF (rho < eps_rho) THEN
         rs = 0.0_dp
      ELSE
         rs = rsfac*rho**(-f13)
      END IF

   END SUBROUTINE calc_rs_single

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param rs ...
! **************************************************************************************************
   SUBROUTINE calc_rs_array(rho, rs)

!   rs parameter : f*rho**(-1/3)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: rs

      INTEGER                                            :: k

      IF (SIZE(rs) < SIZE(rho)) THEN
         CPABORT("Size of array rs too small.")
      END IF

!$OMP PARALLEL DO PRIVATE(k) DEFAULT(NONE) SHARED(rs,eps_rho,rho)
      DO k = 1, SIZE(rs)
         IF (rho(k) < eps_rho) THEN
            rs(k) = 0.0_dp
         ELSE
            rs(k) = rsfac*rho(k)**(-f13)
         END IF
      END DO

   END SUBROUTINE calc_rs_array

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param rs ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE calc_rs_pw(rho, rs, n)

!   rs parameter : f*rho**(-1/3)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: rs
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: k

!$OMP PARALLEL DO PRIVATE(k) SHARED(n,rs,rho,eps_rho) DEFAULT(NONE)
      DO k = 1, n
         IF (rho(k) < eps_rho) THEN
            rs(k) = 0.0_dp
         ELSE
            rs(k) = rsfac*rho(k)**(-f13)
         END IF
      END DO

   END SUBROUTINE calc_rs_pw

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE calc_srs_pw(rho, x, n)

!   rs parameter : f*rho**(-1/3)
!   x = sqrt(rs)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: x
      INTEGER, INTENT(in)                                :: n

      INTEGER                                            :: ip

      CALL calc_rs_pw(rho, x, n)

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE) SHARED(x,n)
      DO ip = 1, n
         x(ip) = SQRT(x(ip))
      END DO

   END SUBROUTINE calc_srs_pw

! **************************************************************************************************
!> \brief ...
!> \param tag ...
!> \param rho ...
!> \param grho ...
!> \param s ...
! **************************************************************************************************
   SUBROUTINE calc_wave_vector(tag, rho, grho, s)

!   wave vector s = |nabla rho| / (2(3pi^2)^1/3 * rho^4/3)

      CHARACTER(len=*), INTENT(IN)                       :: tag
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, grho
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: s

      INTEGER                                            :: ip, n
      REAL(KIND=dp)                                      :: fac

!   TAGS: U: total density, spin wave vector
!         R: spin density, total density wave vector

      fac = 1.0_dp/(2.0_dp*(3.0_dp*pi*pi)**f13)
      IF (tag(1:1) == "u" .OR. tag(1:1) == "U") fac = fac*(2.0_dp)**f13
      IF (tag(1:1) == "r" .OR. tag(1:1) == "R") fac = fac*(2.0_dp)**f13

      n = SIZE(s) !FM it was sizederiv_rho
      !FM IF ( n > SIZE(s) ) &
      !FM   CPABORT("Incompatible array sizes.")
      !FM IF ( n > SIZE(grho) ) &
      !FM   CPABORT("Incompatible array sizes.")

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE) SHARED(rho,eps_rho,s,fac,grho,n)
      DO ip = 1, n
         IF (rho(ip) < eps_rho) THEN
            s(ip) = 0.0_dp
         ELSE
            s(ip) = fac*grho(ip)*rho(ip)**(-f43)
         END IF
      END DO

   END SUBROUTINE calc_wave_vector

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param rhoa ...
!> \param rhob ...
!> \param fx ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE calc_fx_array(n, rhoa, rhob, fx, m)

!   spin interpolation function and derivatives
!
!   f(x) = ( (1+x)^(4/3) + (1-x)^(4/3) - 2 ) / (2^(4/3)-2)
!   df(x) = (4/3)( (1+x)^(1/3) - (1-x)^(1/3) ) / (2^(4/3)-2)
!   d2f(x) = (4/9)( (1+x)^(-2/3) + (1-x)^(-2/3) ) / (2^(4/3)-2)
!   d3f(x) = (-8/27)( (1+x)^(-5/3) - (1-x)^(-5/3) ) / (2^(4/3)-2)
!

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, rhob
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: fx
      INTEGER, INTENT(IN)                                :: m

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: rhoab, x

! number of points
! order of derivatives
!   *** Parameters ***

      IF (m > 3) CPABORT("Order too high.")
!!    IF (.NOT.ASSOCIATED(fx)) THEN
!!       ALLOCATE(fx(n,m+1), STAT=ierr)
!!       IF (ierr /= 0) CALL stop_memory(routineP, "fx", n*m)
!!    ELSE
      IF (SIZE(fx, 1) < n) CPABORT("SIZE(fx,1) too small")
      IF (SIZE(fx, 2) < m) CPABORT("SIZE(fx,2) too small")
!!    END IF

!$OMP PARALLEL DO PRIVATE(ip,x,rhoab) DEFAULT(NONE) SHARED(fx,m,eps_rho,n)
      DO ip = 1, n
         rhoab = rhoa(ip) + rhob(ip)
         IF (rhoab < eps_rho) THEN
            fx(ip, 1:m) = 0.0_dp
         ELSE
            x = (rhoa(ip) - rhob(ip))/rhoab
            IF (x < -1.0_dp) THEN
               IF (m >= 0) fx(ip, 1) = 1.0_dp
               IF (m >= 1) fx(ip, 2) = -f43*fxfac*2.0_dp**f13
               IF (m >= 2) fx(ip, 3) = f13*f43*fxfac/2.0_dp**f23
               IF (m >= 3) fx(ip, 4) = f23*f13*f43*fxfac/2.0_dp**f53
            ELSE IF (x > 1.0_dp) THEN
               IF (m >= 0) fx(ip, 1) = 1.0_dp
               IF (m >= 1) fx(ip, 2) = f43*fxfac*2.0_dp**f13
               IF (m >= 2) fx(ip, 3) = f13*f43*fxfac/2.0_dp**f23
               IF (m >= 3) fx(ip, 4) = -f23*f13*f43*fxfac/2.0_dp**f53
            ELSE
               IF (m >= 0) &
                  fx(ip, 1) = ((1.0_dp + x)**f43 + (1.0_dp - x)**f43 - 2.0_dp)*fxfac
               IF (m >= 1) &
                  fx(ip, 2) = ((1.0_dp + x)**f13 - (1.0_dp - x)**f13)*fxfac*f43
               IF (m >= 2) &
                  fx(ip, 3) = ((1.0_dp + x)**(-f23) + (1.0_dp - x)**(-f23))* &
                              fxfac*f43*f13
               IF (m >= 3) &
                  fx(ip, 4) = ((1.0_dp + x)**(-f53) - (1.0_dp - x)**(-f53))* &
                              fxfac*f43*f13*(-f23)
            END IF
         END IF
      END DO

   END SUBROUTINE calc_fx_array

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param fx ...
!> \param m ...
! **************************************************************************************************
   PURE SUBROUTINE calc_fx_single(rhoa, rhob, fx, m)

!   spin interpolation function and derivatives
!
!   f(x) = ( (1+x)^(4/3) + (1-x)^(4/3) - 2 ) / (2^(4/3)-2)
!   df(x) = (4/3)( (1+x)^(1/3) - (1-x)^(1/3) ) / (2^(4/3)-2)
!   d2f(x) = (4/9)( (1+x)^(-2/3) + (1-x)^(-2/3) ) / (2^(4/3)-2)
!   d3f(x) = (-8/27)( (1+x)^(-5/3) - (1-x)^(-5/3) ) / (2^(4/3)-2)
!

      REAL(KIND=dp), INTENT(IN)                          :: rhoa, rhob
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: fx
      INTEGER, INTENT(IN)                                :: m

      REAL(KIND=dp)                                      :: rhoab, x

      rhoab = rhoa + rhob
      IF (rhoab < eps_rho) THEN
         fx(1:m) = 0.0_dp
      ELSE
         x = (rhoa - rhob)/rhoab
         IF (x < -1.0_dp) THEN
            IF (m >= 0) fx(1) = 1.0_dp
            IF (m >= 1) fx(2) = -f43*fxfac*2.0_dp**f13
            IF (m >= 2) fx(3) = f13*f43*fxfac/2.0_dp**f23
            IF (m >= 3) fx(4) = f23*f13*f43*fxfac/2.0_dp**f53
         ELSE IF (x > 1.0_dp) THEN
            IF (m >= 0) fx(1) = 1.0_dp
            IF (m >= 1) fx(2) = f43*fxfac*2.0_dp**f13
            IF (m >= 2) fx(3) = f13*f43*fxfac/2.0_dp**f23
            IF (m >= 3) fx(4) = -f23*f13*f43*fxfac/2.0_dp**f53
         ELSE
            IF (m >= 0) &
               fx(1) = ((1.0_dp + x)**f43 + (1.0_dp - x)**f43 - 2.0_dp)*fxfac
            IF (m >= 1) &
               fx(2) = ((1.0_dp + x)**f13 - (1.0_dp - x)**f13)*fxfac*f43
            IF (m >= 2) &
               fx(3) = ((1.0_dp + x)**(-f23) + (1.0_dp - x)**(-f23))* &
                       fxfac*f43*f13
            IF (m >= 3) &
               fx(4) = ((1.0_dp + x)**(-f53) - (1.0_dp - x)**(-f53))* &
                       fxfac*f43*f13*(-f23)
         END IF
      END IF

   END SUBROUTINE calc_fx_single

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param z ...
!> \param order ...
! **************************************************************************************************
   PURE SUBROUTINE calc_z(a, b, z, order)

      REAL(KIND=dp), INTENT(IN)                          :: a, b
      REAL(KIND=dp), DIMENSION(0:, 0:), INTENT(OUT)      :: z
      INTEGER, INTENT(IN)                                :: order

      REAL(KIND=dp)                                      :: c, d

      c = a + b

      z(0, 0) = (a - b)/c
      IF (order >= 1) THEN
         d = c*c
         z(1, 0) = 2.0_dp*b/d
         z(0, 1) = -2.0_dp*a/d
      END IF
      IF (order >= 2) THEN
         d = d*c
         z(2, 0) = -4.0_dp*b/d
         z(1, 1) = 2.0_dp*(a - b)/d
         z(0, 2) = 4.0_dp*a/d
      END IF
      IF (order >= 3) THEN
         d = d*c
         z(3, 0) = 12.0_dp*b/d
         z(2, 1) = -4.0_dp*(a - 2.0_dp*b)/d
         z(1, 2) = -4.0_dp*(2.0_dp*a - b)/d
         z(0, 3) = -12.0_dp*a/d
      END IF

   END SUBROUTINE calc_z

END MODULE xc_functionals_utilities

