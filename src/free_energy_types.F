!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief defines types for metadynamics calculation
!> \par History
!>      01.2007 created [tlaino] Teodoro Laino
! **************************************************************************************************
MODULE free_energy_types

   USE input_constants,                 ONLY: do_fe_ac,&
                                              do_fe_ui
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'free_energy_types'

   PUBLIC :: free_energy_type, &
             fe_env_release, &
             fe_env_create, &
             ui_var_type

! **************************************************************************************************
   TYPE ui_var_type
      REAL(KIND=dp), DIMENSION(:), POINTER           :: ss => NULL()
      INTEGER                                        :: icolvar = -1
   END TYPE ui_var_type

! **************************************************************************************************
   TYPE ui_conv_type
      ! Specifying convergence parameters
      INTEGER                                        :: cg_width = -1, max_cg_width = -1
      INTEGER                                        :: cg_points = -1
      REAL(KIND=dp)                                  :: eps_conv = 0.0_dp
      REAL(KIND=dp)                                  :: k_conf_lm = 0.0_dp
      REAL(KIND=dp)                                  :: sw_conf_lm = 0.0_dp
      REAL(KIND=dp)                                  :: vn_conf_lm = 0.0_dp
      LOGICAL                                        :: test_k = .FALSE., &
                                                        test_sw = .FALSE., &
                                                        test_vn = .FALSE.
   END TYPE ui_conv_type

! **************************************************************************************************
   TYPE statistical_type
      ! Collecting coarse grained data
      REAL(KIND=dp), DIMENSION(:), POINTER           :: avg => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER         :: var => NULL()
   END TYPE statistical_type

! **************************************************************************************************
   TYPE free_energy_type
      INTEGER                                        :: ncolvar = -1
      INTEGER                                        :: TYPE = -1
      INTEGER                                        :: nr_points = -1, &
                                                        nr_rejected = -1
      TYPE(ui_conv_type), POINTER                    :: conv_par => NULL()
      TYPE(ui_var_type), POINTER, DIMENSION(:)       :: uivar => NULL()
      TYPE(statistical_type), DIMENSION(:), POINTER  :: cg_data => NULL()
      ! Old data
      REAL(KIND=dp)                                  :: eps_conv = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), POINTER         :: covmx => NULL()

      CHARACTER(len=default_string_length)           :: plumed_input_file = ""
   END TYPE free_energy_type

CONTAINS

! **************************************************************************************************
!> \brief creates the fe_env
!> \param fe_env ...
!> \param fe_section ...
!> \author Teodoro Laino 01.2007
! **************************************************************************************************
   SUBROUTINE fe_env_create(fe_env, fe_section)
      TYPE(free_energy_type), POINTER                    :: fe_env
      TYPE(section_vals_type), POINTER                   :: fe_section

      INTEGER                                            :: i, id_method
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: ui_section, ui_var_section

      CPASSERT(.NOT. ASSOCIATED(fe_env))

      CALL section_vals_get(fe_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(fe_section, "METHOD", i_val=id_method)
         SELECT CASE (id_method)
         CASE (do_fe_ui)
            ALLOCATE (fe_env)
            NULLIFY (fe_env%covmx, fe_env%uivar, fe_env%conv_par, fe_env%cg_data)
            fe_env%type = id_method
            fe_env%nr_points = 0
            fe_env%nr_rejected = 0
            NULLIFY (fe_env%cg_data)
            ui_section => section_vals_get_subs_vals(fe_section, "UMBRELLA_INTEGRATION")
            ui_var_section => section_vals_get_subs_vals(ui_section, "UVAR")
            CALL section_vals_get(ui_var_section, n_repetition=fe_env%ncolvar)
            ! Convergence controlling parameters
            ALLOCATE (fe_env%conv_par)
            fe_env%conv_par%test_k = .FALSE.
            fe_env%conv_par%test_sw = .FALSE.
            fe_env%conv_par%test_vn = .FALSE.
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%COARSE_GRAINED_WIDTH", &
                                      i_val=fe_env%conv_par%cg_width)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%MAX_COARSE_GRAINED_WIDTH", &
                                      i_val=fe_env%conv_par%max_cg_width)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%COARSE_GRAINED_POINTS", &
                                      i_val=fe_env%conv_par%cg_points)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%EPS_CONV", &
                                      r_val=fe_env%conv_par%eps_conv)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%K_CONFIDENCE_LIMIT", &
                                      r_val=fe_env%conv_par%k_conf_lm)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%SW_CONFIDENCE_LIMIT", &
                                      r_val=fe_env%conv_par%sw_conf_lm)
            CALL section_vals_val_get(ui_section, "CONVERGENCE_CONTROL%VN_CONFIDENCE_LIMIT", &
                                      r_val=fe_env%conv_par%vn_conf_lm)
            ! Umbrella Integration variables
            ALLOCATE (fe_env%uivar(fe_env%ncolvar))
            DO i = 1, fe_env%ncolvar
               ! Read Umbrella Integration Variable definition
               CALL section_vals_val_get(ui_var_section, "COLVAR", &
                                         i_val=fe_env%uivar(i)%icolvar, i_rep_section=i)
               NULLIFY (fe_env%uivar(i)%ss)
            END DO
         CASE (do_fe_ac)
            ALLOCATE (fe_env)
            NULLIFY (fe_env%covmx, fe_env%uivar, fe_env%conv_par, fe_env%cg_data)
            ALLOCATE (fe_env%covmx(3, 0))
            fe_env%type = id_method
            CALL section_vals_val_get(fe_section, "ALCHEMICAL_CHANGE%EPS_CONV", r_val=fe_env%eps_conv)
         CASE DEFAULT
            ! Do Nothing
         END SELECT
      END IF
   END SUBROUTINE fe_env_create

! **************************************************************************************************
!> \brief releases the fe_env
!> \param fe_env ...
!> \author Laino Teodoro 01.2007
! **************************************************************************************************
   SUBROUTINE fe_env_release(fe_env)
      TYPE(free_energy_type), POINTER                    :: fe_env

      INTEGER                                            :: i

      IF (ASSOCIATED(fe_env)) THEN
         IF (ASSOCIATED(fe_env%covmx)) THEN
            DEALLOCATE (fe_env%covmx)
         END IF
         IF (ASSOCIATED(fe_env%cg_data)) THEN
            DO i = 1, SIZE(fe_env%cg_data)
               IF (ASSOCIATED(fe_env%cg_data(i)%avg)) THEN
                  DEALLOCATE (fe_env%cg_data(i)%avg)
               END IF
               IF (ASSOCIATED(fe_env%cg_data(i)%var)) THEN
                  DEALLOCATE (fe_env%cg_data(i)%var)
               END IF
            END DO
            DEALLOCATE (fe_env%cg_data)
         END IF
         IF (ASSOCIATED(fe_env%conv_par)) THEN
            DEALLOCATE (fe_env%conv_par)
         END IF
         IF (ASSOCIATED(fe_env%uivar)) THEN
            DO i = 1, SIZE(fe_env%uivar)
               IF (ASSOCIATED(fe_env%uivar(i)%ss)) THEN
                  DEALLOCATE (fe_env%uivar(i)%ss)
               END IF
            END DO
            DEALLOCATE (fe_env%uivar)
         END IF
         DEALLOCATE (fe_env)
      END IF
   END SUBROUTINE fe_env_release

END MODULE free_energy_types
