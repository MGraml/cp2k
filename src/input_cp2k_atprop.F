!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input section for atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! **************************************************************************************************
MODULE input_cp2k_atprop
   USE bibliography,                    ONLY: Kikuchi2009
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_atprop'

   PUBLIC :: create_atprop_section

CONTAINS

! **************************************************************************************************
!> \brief Creates the ATOMIC section
!> \param section the section to create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_atprop_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ATOMIC", &
                          description="Controls the calculation of atomic properties. "// &
                          "Printing is controlled by FORCE_EVAL / PRINT / PROGRAM_RUN_INFO", &
                          repeats=.FALSE., &
                          citations=(/Kikuchi2009/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY", &
                          description="Calculate atomic energies ", &
                          usage="ENERGY {logical}", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_atprop_section

END MODULE input_cp2k_atprop
