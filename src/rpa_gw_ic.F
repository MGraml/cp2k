!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate image charge corrections
!> \par History
!>      06.2019 Moved from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE rpa_gw_ic
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE dbt_api,                         ONLY: &
        dbt_contract, dbt_copy, dbt_copy_matrix_to_tensor, dbt_create, dbt_destroy, dbt_get_block, &
        dbt_get_info, dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, &
        dbt_iterator_stop, dbt_iterator_type, dbt_nblks_total, dbt_pgrid_create, &
        dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_dims_create,&
                                              mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE qs_tensors_types,                ONLY: create_2c_tensor
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_ic'

   PUBLIC :: calculate_ic_correction, apply_ic_corr

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param mat_SinvVSinv ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param unit_nr ...
!> \param print_ic_values ...
!> \param para_env ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_ic_correction(Eigenval, mat_SinvVSinv, &
                                      t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, gw_corr_lev_tot, &
                                      gw_corr_lev_occ, gw_corr_lev_virt, homo, unit_nr, &
                                      print_ic_values, para_env, &
                                      do_alpha, do_beta)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: mat_SinvVSinv
      TYPE(dbt_type)                                     :: t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, unit_nr
      LOGICAL, INTENT(IN)                                :: print_ic_values
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ic_correction'

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, mo_end, mo_start, n_level_gw, &
                                                            n_level_gw_ref
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2, sizes_RI_split
      INTEGER, DIMENSION(2)                              :: pdims
      LOGICAL                                            :: do_closed_shell, my_do_alpha, my_do_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Delta_Sigma_Neaton
      TYPE(dbt_pgrid_type)                               :: pgrid_2d
      TYPE(dbt_type) :: t_3c_overl_nnP_ic_reflected_ctr, t_SinvVSinv, t_SinvVSinv_tmp

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)

      ALLOCATE (Delta_Sigma_Neaton(gw_corr_lev_tot))
      Delta_Sigma_Neaton = 0.0_dp

      mo_start = homo - gw_corr_lev_occ + 1
      mo_end = homo + gw_corr_lev_virt
      CPASSERT(mo_end - mo_start + 1 == gw_corr_lev_tot)

      ALLOCATE (sizes_RI_split(dbt_nblks_total(t_3c_overl_nnP_ic_reflected, 1)))
      CALL dbt_get_info(t_3c_overl_nnP_ic_reflected, blk_size_1=sizes_RI_split)

      CALL dbt_create(mat_SinvVSinv, t_SinvVSinv_tmp)
      CALL dbt_copy_matrix_to_tensor(mat_SinvVSinv, t_SinvVSinv_tmp)
      pdims = 0
      CALL mp_dims_create(para_env%num_pe, pdims)
      CALL dbt_pgrid_create(para_env, pdims, pgrid_2d)
      CALL create_2c_tensor(t_SinvVSinv, dist_1, dist_2, pgrid_2d, sizes_RI_split, sizes_RI_split, &
                            name="(RI|RI)")
      DEALLOCATE (dist_1, dist_2)
      CALL dbt_pgrid_destroy(pgrid_2d)

      CALL dbt_copy(t_SinvVSinv_tmp, t_SinvVSinv)
      CALL dbt_destroy(t_SinvVSinv_tmp)
      CALL dbt_create(t_3c_overl_nnP_ic_reflected, t_3c_overl_nnP_ic_reflected_ctr)
      CALL dbt_contract(0.5_dp, t_SinvVSinv, t_3c_overl_nnP_ic_reflected, &
                        0.0_dp, t_3c_overl_nnP_ic_reflected_ctr, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2, 3], &
                        map_1=[1], map_2=[2, 3])

      CALL trace_ic_gw(t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected_ctr, Delta_Sigma_Neaton, [mo_start, mo_end], para_env)

      Delta_Sigma_Neaton(gw_corr_lev_occ + 1:) = -Delta_Sigma_Neaton(gw_corr_lev_occ + 1:)

      CALL dbt_destroy(t_SinvVSinv)
      CALL dbt_destroy(t_3c_overl_nnP_ic_reflected_ctr)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '----------------------------------------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies of alpha spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-------------------------------------------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'Single-electron energies of beta spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '------------------------------------------------------------------------'
         END IF

         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T3,A)') 'Reference for the ic: Neaton et al., PRL 97, 216405 (2006)'
         WRITE (unit_nr, *) ' '

         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T14,2A)') 'MO     E_n before ic corr           Delta E_ic', &
            '    E_n after ic corr'

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
               n_level_gw_ref, ' ( ', occ_virt, ')  ', &
               Eigenval(n_level_gw_ref)*evolt, &
               Delta_Sigma_Neaton(n_level_gw)*evolt, &
               (Eigenval(n_level_gw_ref) + Delta_Sigma_Neaton(n_level_gw))*evolt

         END DO

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F57.2)') 'IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                      Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                      Eigenval(homo) - &
                                                                      Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F51.2)') 'Alpha IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                            Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                            Eigenval(homo) - &
                                                                            Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F52.2)') 'Beta IC HOMO-LUMO gap (eV)', (Eigenval(homo + 1) + &
                                                                           Delta_Sigma_Neaton(gw_corr_lev_occ + 1) - &
                                                                           Eigenval(homo) - &
                                                                           Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
         END IF

         IF (print_ic_values) THEN

            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'Horizontal list for copying the image charge corrections for use as input:'
            WRITE (unit_nr, '(*(F7.3))') (Delta_Sigma_Neaton(n_level_gw)*evolt, &
                                          n_level_gw=1, gw_corr_lev_tot)

         END IF

      END IF

      Eigenval(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval(homo - gw_corr_lev_occ + 1: &
                                                                              homo + gw_corr_lev_virt) &
                                                                     + Delta_Sigma_Neaton(1:gw_corr_lev_tot)

      CALL timestop(handle)

   END SUBROUTINE calculate_ic_correction

! **************************************************************************************************
!> \brief ...
!> \param t3c_1 ...
!> \param t3c_2 ...
!> \param Delta_Sigma_Neaton ...
!> \param mo_bounds ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE trace_ic_gw(t3c_1, t3c_2, Delta_Sigma_Neaton, mo_bounds, para_env)
      TYPE(dbt_type), INTENT(INOUT)                      :: t3c_1, t3c_2
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: Delta_Sigma_Neaton
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trace_ic_gw'

      INTEGER                                            :: handle, n, n_end, n_end_block, n_start, &
                                                            n_start_block
      INTEGER, DIMENSION(3)                              :: boff, bsize, ind
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_1, block_2
      REAL(KIND=dp), &
         DIMENSION(mo_bounds(2)-mo_bounds(1)+1)          :: Delta_Sigma_Neaton_prv
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(Delta_Sigma_Neaton_prv) == SIZE(Delta_Sigma_Neaton))
      Delta_Sigma_Neaton_prv = 0.0_dp

!$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:Delta_Sigma_Neaton_prv) &
!$OMP SHARED(t3c_1,t3c_2,mo_bounds) &
!$OMP PRIVATE(iter,ind,bsize,boff,block_1,block_2,found) &
!$OMP PRIVATE(n,n_start_block,n_start,n_end_block,n_end)
      CALL dbt_iterator_start(iter, t3c_1)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind, blk_size=bsize, blk_offset=boff)
         IF (ind(2) /= ind(3)) CYCLE
         CALL dbt_get_block(t3c_1, ind, block_1, found)
         CPASSERT(found)
         CALL dbt_get_block(t3c_2, ind, block_2, found)
         IF (.NOT. found) CYCLE

         IF (boff(3) < mo_bounds(1)) THEN
            n_start_block = mo_bounds(1) - boff(3) + 1
            n_start = 1
         ELSE
            n_start_block = 1
            n_start = boff(3) - mo_bounds(1) + 1
         END IF

         IF (boff(3) + bsize(3) - 1 > mo_bounds(2)) THEN
            n_end_block = mo_bounds(2) - boff(3) + 1
            n_end = mo_bounds(2) - mo_bounds(1) + 1
         ELSE
            n_end_block = bsize(3)
            n_end = boff(3) + bsize(3) - mo_bounds(1)
         END IF

         Delta_Sigma_Neaton_prv(n_start:n_end) = &
            Delta_Sigma_Neaton_prv(n_start:n_end) + &
            (/(DOT_PRODUCT(block_1(:, n, n), &
                           block_2(:, n, n)), &
               n=n_start_block, n_end_block)/)
         DEALLOCATE (block_1, block_2)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      Delta_Sigma_Neaton = Delta_Sigma_Neaton + Delta_Sigma_Neaton_prv
      CALL para_env%sum(Delta_Sigma_Neaton)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param ic_corr_list ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param gw_corr_lev_tot ...
!> \param homo ...
!> \param nmo ...
!> \param unit_nr ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE apply_ic_corr(Eigenval, Eigenval_scf, ic_corr_list, &
                            gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                            homo, nmo, unit_nr, do_alpha, do_beta)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval, Eigenval_scf
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ic_corr_list
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, &
                                                            gw_corr_lev_tot, homo, nmo, unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_ic_corr'

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: do_closed_shell, my_do_alpha, my_do_beta
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)

      ! check the number of input image charge corrected levels
      CPASSERT(SIZE(ic_corr_list) == gw_corr_lev_tot)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-----------------------------------------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of alpha spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '--------------------------------------------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of beta spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-------------------------------------------------------------------------'
         END IF

         WRITE (unit_nr, *) ' '

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
               n_level_gw_ref, ' ( ', occ_virt, ')  ', &
               Eigenval(n_level_gw_ref)*evolt, &
               ic_corr_list(n_level_gw)*evolt, &
               (Eigenval(n_level_gw_ref) + ic_corr_list(n_level_gw))*evolt

         END DO

         WRITE (unit_nr, *) ' '

      END IF

      Eigenval(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval(homo - gw_corr_lev_occ + 1: &
                                                                              homo + gw_corr_lev_virt) &
                                                                     + ic_corr_list(1:gw_corr_lev_tot)

      Eigenval_scf(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval_scf(homo - gw_corr_lev_occ + 1: &
                                                                                      homo + gw_corr_lev_virt) &
                                                                         + ic_corr_list(1:gw_corr_lev_tot)

      IF (unit_nr > 0) THEN

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A,F52.2)') 'G0W0 IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A,F46.2)') 'G0W0 Alpha IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A,F47.2)') 'G0W0 Beta IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         END IF

         WRITE (unit_nr, *) ' '

      END IF

      ! for eigenvalue self-consistent GW, all eigenvalues have to be corrected
      ! 1) the occupied; check if there are occupied MOs not being corrected by the IC model
      IF (gw_corr_lev_occ < homo .AND. gw_corr_lev_occ > 0) THEN

         ! calculate average IC contribution for occupied orbitals
         eigen_diff = 0.0_dp

         DO n_level_gw = 1, gw_corr_lev_occ
            eigen_diff = eigen_diff + ic_corr_list(n_level_gw)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_occ

         ! correct the eigenvalues of the occupied orbitals which have not been corrected by the IC model
         DO n_level_gw = 1, homo - gw_corr_lev_occ
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
            Eigenval_scf(n_level_gw) = Eigenval_scf(n_level_gw) + eigen_diff
         END DO

      END IF

      ! 2) the virtual: check if there are virtual orbitals not being corrected by the IC model
      IF (gw_corr_lev_virt < nmo - homo .AND. gw_corr_lev_virt > 0) THEN

         ! calculate average IC correction for virtual orbitals
         eigen_diff = 0.0_dp
         DO n_level_gw = gw_corr_lev_occ + 1, gw_corr_lev_tot
            eigen_diff = eigen_diff + ic_corr_list(n_level_gw)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_virt

         ! correct the eigenvalues of the virtual orbitals which have not been corrected by the IC model
         DO n_level_gw = homo + gw_corr_lev_virt + 1, nmo
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
            Eigenval_scf(n_level_gw) = Eigenval_scf(n_level_gw) + eigen_diff
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_ic_corr

END MODULE rpa_gw_ic
