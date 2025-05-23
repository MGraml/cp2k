!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types used in the flexible partitioning scheme
!> \par History
!>      04.2006 [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE fp_types
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fp_type
   PUBLIC :: fp_env_create, fp_env_release, fp_env_read, fp_env_write

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fp_types'

! **************************************************************************************************
   TYPE fp_type

      ! input related objects
      LOGICAL                        :: use_fp = .FALSE.

      INTEGER                        :: central_atom = -1
      INTEGER, DIMENSION(:), POINTER :: inner_atoms => NULL(), outer_atoms => NULL()
      REAL(KIND=dp)                  :: inner_radius = 0.0_dp, outer_radius = 0.0_dp
      REAL(KIND=dp)                  :: strength = 0.0_dp, smooth_width = 0.0_dp
      LOGICAL                        :: bias = .FALSE.
      REAL(KIND=dp)                  :: temperature = 0.0_dp
      TYPE(section_vals_type), POINTER   :: print_section => NULL()

      ! computed during runs
      INTEGER                        :: i1 = -1, i2 = -1, o1 = -1, o2 = -1
      REAL(KIND=dp)                  :: ri1 = 0.0_dp, ri2 = 0.0_dp, ro1 = 0.0_dp, ro2 = 0.0_dp
      REAL(KIND=dp)                  :: weight = 0.0_dp, comb_weight = 0.0_dp, bias_weight = 0.0_dp
      REAL(KIND=dp)                  :: energy = 0.0_dp, bias_energy = 0.0_dp, restraint_energy = 0.0_dp
   END TYPE fp_type

CONTAINS

! **************************************************************************************************
!> \brief create retain release the flexible partitioning environment
!> \param fp_env ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! **************************************************************************************************
   PURE SUBROUTINE fp_env_create(fp_env)
      TYPE(fp_type), INTENT(OUT)                         :: fp_env

      fp_env%use_fp = .FALSE.
      NULLIFY (fp_env%inner_atoms)
      NULLIFY (fp_env%outer_atoms)
      NULLIFY (fp_env%print_section)

   END SUBROUTINE fp_env_create

! **************************************************************************************************
!> \brief ...
!> \param fp_env ...
! **************************************************************************************************
   SUBROUTINE fp_env_release(fp_env)
      TYPE(fp_type), INTENT(INOUT)                       :: fp_env

      IF (ASSOCIATED(fp_env%inner_atoms)) DEALLOCATE (fp_env%inner_atoms)
      IF (ASSOCIATED(fp_env%outer_atoms)) DEALLOCATE (fp_env%outer_atoms)
      IF (ASSOCIATED(fp_env%print_section)) CALL section_vals_release(fp_env%print_section)
      fp_env%use_fp = .FALSE.

   END SUBROUTINE fp_env_release

! **************************************************************************************************
!> \brief reads the corresponding input section and stores it in the fp_env
!> \param fp_env ...
!> \param fp_section ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE fp_env_read(fp_env, fp_section)
      TYPE(fp_type), INTENT(INOUT)                       :: fp_env
      TYPE(section_vals_type), POINTER                   :: fp_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fp_env_read'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: tmplist

      CALL timeset(routineN, handle)
      CALL section_vals_get(fp_section, explicit=fp_env%use_fp)
      IF (fp_env%use_fp) THEN
         CALL section_vals_val_get(fp_section, "CENTRAL_ATOM", i_val=fp_env%central_atom)

         CALL section_vals_val_get(fp_section, "INNER_ATOMS", i_vals=tmplist)
         ALLOCATE (fp_env%inner_atoms(SIZE(tmplist, 1)))
         fp_env%inner_atoms = tmplist

         CALL section_vals_val_get(fp_section, "OUTER_ATOMS", i_vals=tmplist)
         ALLOCATE (fp_env%outer_atoms(SIZE(tmplist, 1)))
         fp_env%outer_atoms = tmplist

         CALL section_vals_val_get(fp_section, "INNER_RADIUS", r_val=fp_env%inner_radius)
         CALL section_vals_val_get(fp_section, "OUTER_RADIUS", r_val=fp_env%outer_radius)
         CALL section_vals_val_get(fp_section, "STRENGTH", r_val=fp_env%strength)
         CALL section_vals_val_get(fp_section, "SMOOTH_WIDTH", r_val=fp_env%smooth_width)
         CALL section_vals_val_get(fp_section, "BIAS", l_val=fp_env%bias)
         CALL section_vals_val_get(fp_section, "TEMPERATURE", r_val=fp_env%temperature)

         fp_env%print_section => section_vals_get_subs_vals(fp_section, "WEIGHTS")
         CALL section_vals_retain(fp_env%print_section)
      END IF
      CALL timestop(handle)

   END SUBROUTINE fp_env_read

! **************************************************************************************************
!> \brief writes information concerning the fp_env to the output
!> \param fp_env ...
!> \param fp_section ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE fp_env_write(fp_env, fp_section)
      TYPE(fp_type), INTENT(IN)                          :: fp_env
      TYPE(section_vals_type), POINTER                   :: fp_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fp_env_write'

      INTEGER                                            :: handle, output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: PRINT

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (fp_env%use_fp) THEN
         PRINT => section_vals_get_subs_vals(fp_section, "CONTROL")
         output_unit = cp_print_key_unit_nr(logger, PRINT, "", extension=".Log")
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T79,A)") &
               "FP| Flexible partitioning is ", "ON"
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
               "FP| Central atom ", fp_env%central_atom
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
               "FP| number of inner atoms", SIZE(fp_env%inner_atoms, 1)
            WRITE (UNIT=output_unit, FMT="(1(T2,8I8))") fp_env%inner_atoms
            WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
               "FP| number of outer atoms", SIZE(fp_env%outer_atoms, 1)
            WRITE (UNIT=output_unit, FMT="(1(T2,8I8))") fp_env%outer_atoms
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
               "FP| inner radius [a.u.] ", fp_env%inner_radius
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
               "FP| outer radius [a.u.] ", fp_env%outer_radius
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
               "FP| reflecting restraint strength ", fp_env%strength
            IF (fp_env%bias) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,T79,A)") &
                  "FP| Flexible partitioning bias is ", "ON"
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
                  "FP| bias temperature [kT a.u.]", fp_env%temperature
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
                  "FP| smooth width [a.u.] ", fp_env%smooth_width
            ELSE
               WRITE (UNIT=output_unit, FMT="(T2,A,T78,A)") &
                  "FP| Flexible partitioning bias is", "OFF"
            END IF
         END IF
         CALL cp_print_key_finished_output(output_unit, logger, PRINT, "")
      END IF
      CALL timestop(handle)

   END SUBROUTINE fp_env_write

END MODULE fp_types
