!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains a functional that calculates the energy and its derivatives
!>      for the geometry optimizer
!> \par History
!>      01.2008 - Luca Bellucci and Teodoro Laino - Generalizing for Dimer Method.
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! **************************************************************************************************
MODULE gopt_f_types
   USE cell_opt_types,                  ONLY: cell_opt_env_create,&
                                              cell_opt_env_release,&
                                              cell_opt_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE dimer_types,                     ONLY: dimer_env_create,&
                                              dimer_env_release,&
                                              dimer_env_retain,&
                                              dimer_env_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_get_natom,&
                                              force_env_release,&
                                              force_env_retain,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE gopt_param_types,                ONLY: gopt_param_read,&
                                              gopt_param_type
   USE input_constants,                 ONLY: default_cell_method_id,&
                                              default_dimer_method_id,&
                                              default_minimization_method_id,&
                                              default_shellcore_method_id,&
                                              default_ts_method_id
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_list_types,             ONLY: particle_list_type
   USE space_groups_types,              ONLY: release_spgr_type,&
                                              spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_f_types'

   PUBLIC :: gopt_f_type, gopt_f_create, gopt_f_retain, gopt_f_release

! **************************************************************************************************
!> \brief calculates the potential energy of a system, and its derivatives
!> \par History
!>      none
! **************************************************************************************************
   TYPE gopt_f_type
      INTEGER                                  :: ref_count = 0
      INTEGER                                  :: nfree = 0
 INTEGER                                  :: type_id=default_cell_method_id, ts_method_id=0, cell_method_id=0, shellcore_method_id=0
      LOGICAL                                  :: dimer_rotation = .FALSE., do_line_search = .FALSE., eval_opt_geo = .FALSE.
      CHARACTER(LEN=default_string_length)     :: label = "", tag = ""
      TYPE(force_env_type), POINTER            :: force_env => NULL()
      TYPE(global_environment_type), POINTER   :: globenv => NULL()
      ! Motion section must be references only for IO of the MOTION%PRINT..
      TYPE(section_vals_type), POINTER         :: motion_section => NULL(), geo_section => NULL()
      TYPE(dimer_env_type), POINTER            :: dimer_env => NULL()
      TYPE(gopt_f_type), POINTER               :: gopt_dimer_env => NULL()
      TYPE(gopt_param_type), POINTER           :: gopt_dimer_param => NULL()
      TYPE(cell_opt_env_type), POINTER         :: cell_env => NULL()
      TYPE(spgr_type), POINTER                 :: spgr => NULL()
      REAL(KIND=dp), DIMENSION(3, 3)           :: h_ref = 0.0_dp
      LOGICAL                                  :: require_consistent_energy_force = .FALSE.
   END TYPE gopt_f_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gopt_env the geometry optimization environment to be created
!>      force_env:
!> \param gopt_param ...
!> \param force_env ...
!> \param globenv ...
!> \param geo_opt_section ...
!> \param eval_opt_geo ...
!> \par History
!>      none
! **************************************************************************************************
   RECURSIVE SUBROUTINE gopt_f_create(gopt_env, gopt_param, force_env, globenv, geo_opt_section, &
                                      eval_opt_geo)

      TYPE(gopt_f_type), POINTER                         :: gopt_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_opt_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: eval_opt_geo

      INTEGER                                            :: natom, nshell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles, shell_particles
      TYPE(section_vals_type), POINTER                   :: dimer_section, rot_opt_section

      CPASSERT(.NOT. ASSOCIATED(gopt_env))
      ALLOCATE (gopt_env)
      nshell = 0

      NULLIFY (gopt_env%dimer_env, gopt_env%gopt_dimer_env, gopt_env%gopt_dimer_param, gopt_env%cell_env, gopt_env%spgr)
      gopt_env%ref_count = 1
      gopt_env%dimer_rotation = .FALSE.
      gopt_env%do_line_search = .FALSE.
      ALLOCATE (gopt_env%spgr)
      CALL force_env_retain(force_env)
      gopt_env%force_env => force_env
      gopt_env%motion_section => section_vals_get_subs_vals(force_env%root_section, "MOTION")
      gopt_env%geo_section => geo_opt_section
      gopt_env%globenv => globenv
      gopt_env%eval_opt_geo = .TRUE.
      IF (PRESENT(eval_opt_geo)) gopt_env%eval_opt_geo = eval_opt_geo
      gopt_env%require_consistent_energy_force = .TRUE.

      CALL force_env_get(force_env, subsys=subsys)
      gopt_env%type_id = gopt_param%type_id
      SELECT CASE (gopt_env%type_id)
      CASE (default_ts_method_id, default_minimization_method_id)
         CALL cp_subsys_get(subsys, &
                            particles=particles, &
                            shell_particles=shell_particles)
         IF (ASSOCIATED(shell_particles)) nshell = shell_particles%n_els
         ! The same number of shell and core particles is assumed
         gopt_env%nfree = particles%n_els + nshell
         gopt_env%label = "GEO_OPT"
         gopt_env%tag = "GEOMETRY"
         SELECT CASE (gopt_param%type_id)
         CASE (default_ts_method_id)
            gopt_env%ts_method_id = gopt_param%ts_method_id
            SELECT CASE (gopt_param%ts_method_id)
            CASE (default_dimer_method_id)
               ! For the Dimer method we use the same framework of geometry optimizers
               ! already defined for cp2k..
               natom = force_env_get_natom(force_env)
               dimer_section => section_vals_get_subs_vals(geo_opt_section, "TRANSITION_STATE%DIMER")
               CALL dimer_env_create(gopt_env%dimer_env, subsys, globenv, dimer_section)

               ! Setup the GEO_OPT environment for the rotation of the Dimer
               rot_opt_section => section_vals_get_subs_vals(dimer_section, "ROT_OPT")
               ALLOCATE (gopt_env%gopt_dimer_param)
               CALL gopt_param_read(gopt_env%gopt_dimer_param, rot_opt_section, &
                                    type_id=default_minimization_method_id)
               gopt_env%gopt_dimer_param%type_id = default_ts_method_id

               CALL gopt_f_create(gopt_env%gopt_dimer_env, gopt_env%gopt_dimer_param, force_env=force_env, &
                                  globenv=globenv, geo_opt_section=rot_opt_section, eval_opt_geo=eval_opt_geo)
               CALL dimer_env_retain(gopt_env%dimer_env)
               gopt_env%gopt_dimer_env%dimer_env => gopt_env%dimer_env
               gopt_env%gopt_dimer_env%label = "ROT_OPT"
               gopt_env%gopt_dimer_env%dimer_rotation = .TRUE.
            END SELECT
         END SELECT
      CASE (default_cell_method_id)
         gopt_env%nfree = 6
         gopt_env%label = "CELL_OPT"
         gopt_env%tag = "  CELL  "
         gopt_env%cell_method_id = gopt_param%cell_method_id
         ALLOCATE (gopt_env%cell_env)
         CALL cell_opt_env_create(gopt_env%cell_env, force_env, gopt_env%geo_section)
      CASE (default_shellcore_method_id)
         gopt_env%nfree = subsys%shell_particles%n_els
         gopt_env%label = "SHELL_OPT"
         gopt_env%tag = "  SHELL-CORE  "
         gopt_env%shellcore_method_id = gopt_param%shellcore_method_id
      END SELECT
   END SUBROUTINE gopt_f_create

! **************************************************************************************************
!> \brief ...
!> \param gopt_env the geometry optimization environment to retain
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE gopt_f_retain(gopt_env)
      TYPE(gopt_f_type), POINTER                         :: gopt_env

      CPASSERT(ASSOCIATED(gopt_env))
      CPASSERT(gopt_env%ref_count > 0)
      gopt_env%ref_count = gopt_env%ref_count + 1
   END SUBROUTINE gopt_f_retain

! **************************************************************************************************
!> \brief ...
!> \param gopt_env the geometry optimization environment to release
!> \par History
!>      none
! **************************************************************************************************
   RECURSIVE SUBROUTINE gopt_f_release(gopt_env)
      TYPE(gopt_f_type), POINTER                         :: gopt_env

      IF (ASSOCIATED(gopt_env)) THEN
         CPASSERT(gopt_env%ref_count > 0)
         gopt_env%ref_count = gopt_env%ref_count - 1
         IF (gopt_env%ref_count == 0) THEN
            CALL force_env_release(gopt_env%force_env)
            NULLIFY (gopt_env%force_env, &
                     gopt_env%globenv, &
                     gopt_env%motion_section, &
                     gopt_env%geo_section)
            IF (ASSOCIATED(gopt_env%cell_env)) THEN
               CALL cell_opt_env_release(gopt_env%cell_env)
               DEALLOCATE (gopt_env%cell_env)
            END IF
            CALL dimer_env_release(gopt_env%dimer_env)
            CALL gopt_f_release(gopt_env%gopt_dimer_env)
            IF (ASSOCIATED(gopt_env%gopt_dimer_param)) DEALLOCATE (gopt_env%gopt_dimer_param)
            CALL release_spgr_type(gopt_env%spgr)
            DEALLOCATE (gopt_env)
         END IF
      END IF
   END SUBROUTINE gopt_f_release

END MODULE gopt_f_types
