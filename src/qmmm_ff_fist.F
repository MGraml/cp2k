!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author teo
! **************************************************************************************************
MODULE qmmm_ff_fist

   USE kinds,                           ONLY: default_string_length
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_ff_fist'

   PUBLIC :: qmmm_ff_precond_only_qm

CONTAINS

! **************************************************************************************************
!> \brief This function handles the atom names and modifies
!>      the "_QM_" prefix, in order to find the parameters
!>      and return .TRUE. if all input labels belong to QM atoms
!> \param id1 ...
!> \param id2 ...
!> \param id3 ...
!> \param id4 ...
!> \param is_link ...
!> \return ...
!> \par History
!>      11.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   FUNCTION qmmm_ff_precond_only_qm(id1, id2, id3, id4, is_link) RESULT(only_qm)
      CHARACTER(LEN=default_string_length), &
         INTENT(INOUT)                                   :: id1
      CHARACTER(LEN=default_string_length), &
         INTENT(INOUT), OPTIONAL                         :: id2, id3, id4
      LOGICAL, INTENT(OUT), OPTIONAL                     :: is_link
      LOGICAL                                            :: only_qm

      CHARACTER(LEN=default_string_length)               :: tmp
      INTEGER                                            :: j, my_index
      LOGICAL                                            :: my_link

      only_qm = .FALSE.
      my_link = .FALSE.
      DO WHILE (INDEX(id1, "_QM_") /= 0)
         my_link = qmmm_ff_precond_only_link(id1) .OR. my_link
         my_index = INDEX(id1, "_QM_") + LEN_TRIM("_QM_")
         only_qm = .TRUE.
         tmp = TRIM(id1(my_index:))
         clean_string_1: DO j = 1, default_string_length
            id1(j:j) = " "
         END DO clean_string_1
         id1 = TRIM(tmp)
      END DO

      IF (PRESENT(id2)) THEN
         IF (INDEX(id2, "_QM_") == 0) only_qm = .FALSE.
         DO WHILE (INDEX(id2, "_QM_") /= 0)
            my_link = qmmm_ff_precond_only_link(id2) .OR. my_link
            my_index = INDEX(id2, "_QM_") + LEN_TRIM("_QM_")
            tmp = TRIM(id2(my_index:))
            clean_string_2: DO j = 1, default_string_length
               id2(j:j) = " "
            END DO clean_string_2
            id2 = TRIM(tmp)
         END DO
      END IF

      IF (PRESENT(id3)) THEN
         IF (INDEX(id3, "_QM_") == 0) only_qm = .FALSE.
         DO WHILE (INDEX(id3, "_QM_") /= 0)
            my_link = qmmm_ff_precond_only_link(id3) .OR. my_link
            my_index = INDEX(id3, "_QM_") + LEN_TRIM("_QM_")
            tmp = TRIM(id3(my_index:))
            clean_string_3: DO j = 1, default_string_length
               id3(j:j) = " "
            END DO clean_string_3
            id3 = TRIM(tmp)
         END DO
      END IF

      IF (PRESENT(id4)) THEN
         IF (INDEX(id4, "_QM_") == 0) only_qm = .FALSE.
         DO WHILE (INDEX(id4, "_QM_") /= 0)
            my_link = qmmm_ff_precond_only_link(id4) .OR. my_link
            my_index = INDEX(id4, "_QM_") + LEN_TRIM("_QM_")
            tmp = TRIM(id4(my_index:))
            clean_string_4: DO j = 1, default_string_length
               id4(j:j) = " "
            END DO clean_string_4
            id4 = TRIM(tmp)
         END DO
      END IF

      IF (PRESENT(is_link)) is_link = my_link

   END FUNCTION qmmm_ff_precond_only_qm

! **************************************************************************************************
!> \brief ...
!> \param id1 ...
!> \return ...
! **************************************************************************************************
   FUNCTION qmmm_ff_precond_only_link(id1) RESULT(is_link)
      CHARACTER(LEN=default_string_length), &
         INTENT(INOUT)                                   :: id1
      LOGICAL                                            :: is_link

      CHARACTER(LEN=default_string_length)               :: tmp
      INTEGER                                            :: j, my_index

      is_link = .FALSE.
      DO WHILE (INDEX(id1, "_LNK") /= 0)
         my_index = INDEX(id1, "_LNK") + 1
         my_index = INDEX(id1(my_index:), "_QM_") + my_index - 1
         is_link = .TRUE.
         tmp = TRIM(id1(my_index:))
         clean_string_1: DO j = 1, default_string_length
            id1(j:j) = " "
         END DO clean_string_1
         id1 = TRIM(tmp)
      END DO

   END FUNCTION qmmm_ff_precond_only_link

END MODULE qmmm_ff_fist
