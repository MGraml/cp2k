!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RT-BSE specific RI-RS kernels (V_grid, W0_grid) and contraction routines.
!>        Independent of which GW variant produced bs_env%fm_W_MIC_freq_zero: this module
!>        builds the RI-RS grid on demand (mat_phi_mu_l, mat_Z_lP) if missing, and assembles
!>        the V_grid and W0_grid kernels that the linearized RT-BSE propagation contracts
!>        against rho/S^-1 to obtain Sigma and the Hartree potential.
!> \author Maximilian Graml
!> \date 05.2026
! **************************************************************************************************
MODULE rt_bse_ri_rs
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_release, dbcsr_scale, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE gw_large_cell_Gamma_ri_rs,       ONLY: atomic_basis_at_grid_point,&
                                              compute_coeff_Z_lP,&
                                              contract_A_B_A,&
                                              hadamard_product_inplace,&
                                              release_dbcsr_topology_and_matrices,&
                                              ri_rs_grid_assembler,&
                                              setup_square_topology
   USE gw_large_cell_gamma,             ONLY: multiply_fm_W_MIC_time_with_Minv_Gamma
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_ri_rs'

   PUBLIC :: rt_bse_ri_rs_ensure_grid, &
             rt_bse_ri_rs_ensure_V_grid, &
             rt_bse_ri_rs_ensure_W0_grid, &
             rt_bse_ri_rs_ensure_psi, &
             compute_sigma_ri_rs, &
             compute_sigma_ri_rs_complex, &
             compute_sigma_ri_rs_mo_complex, &
             compute_hartree_ri_rs, &
             compute_hartree_ri_rs_complex, &
             compute_hartree_ri_rs_mo, &
             compute_hartree_ri_rs_mo_shared_complex

CONTAINS

! **************************************************************************************************
!> \brief Make sure mat_phi_mu_l and mat_Z_lP are populated in memory.
!>        If GW was run with RI-RS the grid is already built; otherwise build it here so the
!>        AO-RI GW + RI-RS RT-BSE combination is possible.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_grid'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      t1 = m_walltime()

      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)
      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      bs_env%ri_rs%grid_built = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Built RI-RS grid for RT-BSE (no GW_RI_RS used),', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_grid

! **************************************************************************************************
!> \brief Build V_aux(P,Q) (M^-1-sandwiched truncated Coulomb in the RI basis); the grid kernel
!>        Z V Z^T is applied factorized per call. Needed by RT-BSE RI-RS Hartree.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_V_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_V_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%V_grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL rt_bse_ri_rs_ensure_grid(bs_env, qs_env)
      t1 = m_walltime()

      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      ! Apply M^-1 sandwich to match the W^MIC convention; same scale as W^c when both used.
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      ! Store the M^-1-sandwiched RI-basis Coulomb; the grid kernel Z V Z^T is applied factorized.
      CALL dbcsr_create(bs_env%ri_rs%mat_V_aux_rtbse, "V_aux_rtbse", dist_aux_aux, &
                        dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), bs_env%ri_rs%mat_V_aux_rtbse, &
                            keep_sparsity=.FALSE.)

      bs_env%ri_rs%V_grid_built = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Precomputed RT-BSE RI-RS V_aux kernel,', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_V_grid

! **************************************************************************************************
!> \brief Build W0_grid(l,l') = sum_PQ Z_lP (V+W^c)_PQ Z_l'Q (full W at omega=0 on the grid).
!>        Needed by RT-BSE RI-RS SEX/COH. Reuses bs_env%fm_W_MIC_freq_zero which must already
!>        contain M^-1 W^c(w=0) M^-1 (built by either GW path under the BSE rtp_method gate).
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_W0_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_W0_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_V_aux, matrix_W_aux

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%W0_grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL rt_bse_ri_rs_ensure_grid(bs_env, qs_env)
      t1 = m_walltime()

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)
      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux_rtbse_W0", dist_aux_aux, dbcsr_type_no_symmetry, &
                        blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      CALL dbcsr_create(matrix_W_aux, "W_aux_rtbse", dist_aux_aux, dbcsr_type_no_symmetry, &
                        blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, matrix_W_aux, keep_sparsity=.FALSE.)
      CALL dbcsr_add(matrix_W_aux, matrix_V_aux, 1.0_dp, 1.0_dp)

      CALL dbcsr_create(bs_env%ri_rs%mat_W0_grid_rtbse, "W0_grid_rtbse", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_Z_lP, matrix_W_aux, &
                          bs_env%ri_rs%mat_W0_grid_rtbse, bs_env%eps_filter)

      bs_env%ri_rs%W0_grid_built = .TRUE.
      bs_env%ri_rs%rtbse_kernels_ready = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Precomputed RT-BSE RI-RS W0_grid kernel,', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux, &
                                               m1=matrix_V_aux, m2=matrix_W_aux)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_W0_grid

! **************************************************************************************************
!> \brief Build per-spin MO-active collocation psi_sigma = phi * C_active(sigma) (grid x mo_active).
!>        Idempotent via psi_built; C_active is propagation-constant, one build per run.
!> \param bs_env ...
!> \param C_active per-spin truncated MO coefficient slabs (n_ao x mo_active)
!> \param n_spin ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_psi(bs_env, C_active, n_spin)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: C_active
      INTEGER, INTENT(IN)                                :: n_spin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_psi'

      INTEGER                                            :: chunk, handle, i_blk, i_spin, mo_active, &
                                                            n_mo_chunks, npcol, nprow
      INTEGER, DIMENSION(:), POINTER                     :: ao_blk, col_dist_ao, col_dist_mo, &
                                                            grid_blk, mo_blk_sizes, &
                                                            row_dist_ao_mapped, row_dist_grid
      TYPE(dbcsr_distribution_type)                      :: dist_C, dist_phi, dist_psi
      TYPE(dbcsr_type)                                   :: matrix_C

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%psi_built) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CPASSERT(bs_env%ri_rs%grid_built)

      CALL cp_fm_get_info(C_active(1), ncol_global=mo_active)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, distribution=dist_phi, &
                          row_blk_size=grid_blk, col_blk_size=ao_blk)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, col_dist=col_dist_ao)
      nprow = MAXVAL(row_dist_grid) + 1
      npcol = MAXVAL(col_dist_ao) + 1

      ! mo column chunks (chunk_size_dbcsr convention) and AO row mapping (same trick as
      ! setup_square_topology 'COL')
      chunk = bs_env%ri_rs%chunk_size_dbcsr
      n_mo_chunks = (mo_active + chunk - 1)/chunk
      ALLOCATE (mo_blk_sizes(n_mo_chunks), col_dist_mo(n_mo_chunks))
      DO i_blk = 1, n_mo_chunks
         mo_blk_sizes(i_blk) = MIN(chunk, mo_active - (i_blk - 1)*chunk)
         col_dist_mo(i_blk) = MOD(i_blk - 1, npcol)
      END DO
      ALLOCATE (row_dist_ao_mapped(SIZE(ao_blk)))
      DO i_blk = 1, SIZE(ao_blk)
         row_dist_ao_mapped(i_blk) = MOD(i_blk - 1, nprow)
      END DO

      CALL dbcsr_distribution_new(dist_C, template=dist_phi, row_dist=row_dist_ao_mapped, &
                                  col_dist=col_dist_mo)
      CALL dbcsr_distribution_new(dist_psi, template=dist_phi, row_dist=row_dist_grid, &
                                  col_dist=col_dist_mo)

      ALLOCATE (bs_env%ri_rs%mat_psi_m_l(n_spin))
      DO i_spin = 1, n_spin
         CALL dbcsr_create(matrix_C, "C_active_dbcsr", dist_C, dbcsr_type_no_symmetry, &
                           ao_blk, mo_blk_sizes)
         CALL copy_fm_to_dbcsr(C_active(i_spin), matrix_C, keep_sparsity=.FALSE.)
         CALL dbcsr_create(bs_env%ri_rs%mat_psi_m_l(i_spin), "mat_psi_m_l", dist_psi, &
                           dbcsr_type_no_symmetry, grid_blk, mo_blk_sizes)
         CALL dbcsr_multiply("N", "N", 1.0_dp, bs_env%ri_rs%mat_phi_mu_l, matrix_C, &
                             0.0_dp, bs_env%ri_rs%mat_psi_m_l(i_spin), &
                             filter_eps=bs_env%eps_filter)
         CALL dbcsr_release(matrix_C)
      END DO

      CALL dbcsr_distribution_release(dist_C)
      CALL dbcsr_distribution_release(dist_psi)
      DEALLOCATE (mo_blk_sizes, col_dist_mo, row_dist_ao_mapped)
      bs_env%ri_rs%psi_built = .TRUE.

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_psi

! **************************************************************************************************
!> \brief Real-input self-energy via RI-RS grid kernel:
!>        Sigma_AO(mu,nu) = prefactor * sum_l phi(l,mu) [ rho_grid(l,l') . W0_grid(l,l') ] phi(l',nu)
!>        where rho_grid = phi * rho_AO * phi^T (l x l).
!>        Mirrors the AO-RI get_sigma(rtbse_env, sigma_fm, prefactor, rho_fm) API.
!> \param bs_env ...
!> \param sigma_AO_fm result, AO x AO
!> \param prefactor scaling applied to the final result
!> \param rho_AO_fm input density-like matrix, AO x AO
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs(bs_env, sigma_AO_fm, prefactor, rho_AO_fm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: sigma_AO_fm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_AO_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL sigma_ri_rs_worker(bs_env, sigma_AO_fm, prefactor, rho_AO_fm, &
                              bs_env%ri_rs%mat_phi_mu_l)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs

! **************************************************************************************************
!> \brief SEX kernel body, generic in the collocation matrix (phi: AO domain, psi: MO-active).
!> \param bs_env ...
!> \param sigma_fm result on the colloc column domain (square)
!> \param prefactor scaling applied to the final result
!> \param rho_fm input density-like matrix on the colloc column domain (square)
!> \param mat_colloc collocation matrix (grid x domain)
! **************************************************************************************************
   SUBROUTINE sigma_ri_rs_worker(bs_env, sigma_fm, prefactor, rho_fm, mat_colloc)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: sigma_fm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_fm
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_colloc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sigma_ri_rs_worker'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_grid, blk_in, dist_col_grid, &
                                                            dist_row_in
      TYPE(dbcsr_distribution_type)                      :: dist_grid_grid, dist_in_in
      TYPE(dbcsr_type)                                   :: matrix_rho, matrix_rho_grid, matrix_Sigma

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%W0_grid_built)

      CALL setup_square_topology(mat_colloc, 'COL', dist_in_in, blk_in, dist_row_in)
      CALL setup_square_topology(mat_colloc, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)

      CALL dbcsr_create(matrix_rho, "rho_in_ri_rs", dist_in_in, dbcsr_type_no_symmetry, &
                        blk_in, blk_in)
      CALL copy_fm_to_dbcsr(rho_fm, matrix_rho, keep_sparsity=.FALSE.)

      CALL dbcsr_create(matrix_rho_grid, "rho_grid_ri_rs", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", mat_colloc, matrix_rho, &
                          matrix_rho_grid, bs_env%eps_filter)

      CALL hadamard_product_inplace(matrix_rho_grid, bs_env%ri_rs%mat_W0_grid_rtbse, &
                                    1.0_dp)

      CALL dbcsr_create(matrix_Sigma, template=matrix_rho)
      CALL contract_A_B_A("T", "N", mat_colloc, matrix_rho_grid, &
                          matrix_Sigma, bs_env%eps_filter)

      CALL dbcsr_scale(matrix_Sigma, prefactor)
      CALL copy_dbcsr_to_fm(matrix_Sigma, sigma_fm)

      CALL release_dbcsr_topology_and_matrices(dist=dist_in_in, mapped_dist=dist_row_in, &
                                               m1=matrix_rho, m2=matrix_Sigma)
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=matrix_rho_grid)

      CALL timestop(handle)

   END SUBROUTINE sigma_ri_rs_worker

! **************************************************************************************************
!> \brief Complex-input self-energy via RI-RS: applies compute_sigma_ri_rs to Re and Im parts
!>        independently and reassembles. Required for non-Hermitian Delta-rho inputs in TDA
!>        (kernel must be complex-linear, see rt_bse_linearized notes).
!> \param bs_env ...
!> \param sigma_AO_cfm result, AO x AO (complex)
!> \param prefactor scaling applied to the final result
!> \param rho_AO_cfm input AO x AO complex matrix
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs_complex(bs_env, sigma_AO_cfm, prefactor, rho_AO_cfm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: sigma_AO_cfm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_AO_cfm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_sigma

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_rho, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_sigma, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargetr=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho)
      CALL cp_fm_to_cfm(msourcer=fm_sigma, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargeti=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho)
      CALL cp_fm_to_cfm(msourcei=fm_sigma, mtarget=sigma_AO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), sigma_AO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_sigma)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs_complex

! **************************************************************************************************
!> \brief Hartree potential via RI-RS grid kernel, evaluated diagonal-only (no grid x grid):
!>        V_H_AO(la,si) = sum_l phi(l,la) v_l phi(l,si),  v_l = sum_l' V_grid(l,l') n_l',
!>        n_l = (phi * rho_AO * phi^T)_{ll}, built as the row-dot n_l = sum_n (phi*rho)(l,n) phi(l,n).
!> [INTERNAL] dev_notes/rt_bse_equations.md Eq. VHt_rirs_fact
!>        Real input, real output. Used as building block for the complex variant.
!> \param bs_env ...
!> \param rho_AO_fm input AO x AO density matrix
!> \param V_H_AO_fm output AO x AO Hartree potential
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs(bs_env, rho_AO_fm, V_H_AO_fm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_AO_fm
      TYPE(cp_fm_type), INTENT(INOUT)                    :: V_H_AO_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs'

      INTEGER                                            :: handle, n_grid
      REAL(KIND=dp), ALLOCATABLE                         :: n_vec(:), v_vec(:)

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, nfullrows_total=n_grid)
      ALLOCATE (n_vec(n_grid), v_vec(n_grid))
      n_vec(:) = 0.0_dp

      CALL hartree_density_ri_rs(bs_env, rho_AO_fm, bs_env%ri_rs%mat_phi_mu_l, n_vec, 1.0_dp)
      CALL hartree_potential_ri_rs(bs_env, n_vec, v_vec)
      CALL hartree_backtransform_ri_rs(bs_env, v_vec, V_H_AO_fm, bs_env%ri_rs%mat_phi_mu_l)

      DEALLOCATE (n_vec, v_vec)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs

! **************************************************************************************************
!> \brief Stage 1: accumulate the grid-point density n_vec += fac * diag(colloc*rho*colloc^T).
!> \param bs_env ...
!> \param rho_fm input density on the colloc column domain (square)
!> \param mat_colloc collocation matrix (grid x domain)
!> \param n_vec replicated grid vector, accumulated into (caller zeroes)
!> \param fac accumulation factor (spin degeneracy in the shared cross-spin sum)
! **************************************************************************************************
   SUBROUTINE hartree_density_ri_rs(bs_env, rho_fm, mat_colloc, n_vec, fac)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_fm
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_colloc
      REAL(KIND=dp), INTENT(INOUT)                       :: n_vec(:)
      REAL(KIND=dp), INTENT(IN)                          :: fac

      INTEGER                                            :: n_grid
      INTEGER, DIMENSION(:), POINTER                     :: blk_in, dist_row_in
      REAL(KIND=dp), ALLOCATABLE                         :: tmp_vec(:)
      TYPE(dbcsr_distribution_type)                      :: dist_in_in
      TYPE(dbcsr_type)                                   :: matrix_rho, matrix_T

      CALL setup_square_topology(mat_colloc, 'COL', dist_in_in, blk_in, dist_row_in)
      CALL dbcsr_create(matrix_rho, "rho_hartree_in", dist_in_in, dbcsr_type_no_symmetry, &
                        blk_in, blk_in)
      CALL copy_fm_to_dbcsr(rho_fm, matrix_rho, keep_sparsity=.FALSE.)

      CALL dbcsr_get_info(mat_colloc, nfullrows_total=n_grid)
      ALLOCATE (tmp_vec(n_grid))

      ! diagonal-only density on the grid: T = colloc*rho, n_l = sum_n T(l,n)*colloc(l,n)
      CALL dbcsr_create(matrix_T, template=mat_colloc)
      CALL dbcsr_multiply("N", "N", 1.0_dp, mat_colloc, matrix_rho, &
                          0.0_dp, matrix_T, filter_eps=bs_env%eps_filter)
      CALL grid_row_dot(matrix_T, mat_colloc, tmp_vec, bs_env%para_env)
      CALL dbcsr_release(matrix_T)

      n_vec(:) = n_vec(:) + fac*tmp_vec(:)
      DEALLOCATE (tmp_vec)

      CALL release_dbcsr_topology_and_matrices(dist=dist_in_in, mapped_dist=dist_row_in, &
                                               m1=matrix_rho)

   END SUBROUTINE hartree_density_ri_rs

! **************************************************************************************************
!> \brief Stage 2: v = Z * (V_aux * (Z^T * n)) -- factorized grid Coulomb.
!> \param bs_env ...
!> \param n_vec replicated grid density vector (in)
!> \param v_vec replicated grid potential vector (out)
! **************************************************************************************************
   SUBROUTINE hartree_potential_ri_rs(bs_env, n_vec, v_vec)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: n_vec(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: v_vec(:)

      REAL(KIND=dp), ALLOCATABLE                         :: u_ri(:), w_ri(:)

      ALLOCATE (w_ri(bs_env%n_RI), u_ri(bs_env%n_RI))
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_Z_lP, n_vec, w_ri, bs_env%para_env, &
                                   transposed=.TRUE.)
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_V_aux_rtbse, w_ri, u_ri, bs_env%para_env)
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_Z_lP, u_ri, v_vec, bs_env%para_env)
      DEALLOCATE (w_ri, u_ri)

   END SUBROUTINE hartree_potential_ri_rs

! **************************************************************************************************
!> \brief Stage 3: back-transform V(m,n) = sum_l colloc(l,m) v_l colloc(l,n).
!> \param bs_env ...
!> \param v_vec replicated grid potential vector (in)
!> \param V_fm result on the colloc column domain (square)
!> \param mat_colloc collocation matrix (grid x domain)
! **************************************************************************************************
   SUBROUTINE hartree_backtransform_ri_rs(bs_env, v_vec, V_fm, mat_colloc)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: v_vec(:)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: V_fm
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_colloc

      INTEGER, DIMENSION(:), POINTER                     :: blk_in, dist_row_in
      TYPE(dbcsr_distribution_type)                      :: dist_in_in
      TYPE(dbcsr_type)                                   :: matrix_Phi, matrix_V

      CALL setup_square_topology(mat_colloc, 'COL', dist_in_in, blk_in, dist_row_in)

      ! Phi = rowscale(colloc, v); V = colloc^T * Phi
      CALL dbcsr_create(matrix_Phi, template=mat_colloc)
      CALL dbcsr_copy(matrix_Phi, mat_colloc)
      CALL dbcsr_scale_rows_replicated(matrix_Phi, v_vec)

      CALL dbcsr_create(matrix_V, "V_H_out_ri_rs", dist_in_in, dbcsr_type_no_symmetry, &
                        blk_in, blk_in)
      CALL dbcsr_multiply("T", "N", 1.0_dp, mat_colloc, matrix_Phi, &
                          0.0_dp, matrix_V, filter_eps=bs_env%eps_filter)
      CALL dbcsr_release(matrix_Phi)

      CALL copy_dbcsr_to_fm(matrix_V, V_fm)

      CALL release_dbcsr_topology_and_matrices(dist=dist_in_in, mapped_dist=dist_row_in, &
                                               m1=matrix_V)

   END SUBROUTINE hartree_backtransform_ri_rs

! **************************************************************************************************
!> \brief Complex-input Hartree potential via RI-RS. Re/Im split: feed each part to the real
!>        compute_hartree_ri_rs and reassemble. Real-input Hartree on a non-Hermitian input
!>        would silently drop Im and break Hermitian conjugacy of OV+VO contributions in TDA.
!> \param bs_env ...
!> \param rho_AO_cfm input AO x AO complex density-like matrix
!> \param V_H_AO_cfm output AO x AO complex Hartree potential
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs_complex(bs_env, rho_AO_cfm, V_H_AO_cfm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_AO_cfm
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: V_H_AO_cfm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_v

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_rho, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_v, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargetr=fm_rho)
      CALL compute_hartree_ri_rs(bs_env, fm_rho, fm_v)
      CALL cp_fm_to_cfm(msourcer=fm_v, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargeti=fm_rho)
      CALL compute_hartree_ri_rs(bs_env, fm_rho, fm_v)
      CALL cp_fm_to_cfm(msourcei=fm_v, mtarget=V_H_AO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), V_H_AO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_v)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs_complex

! **************************************************************************************************
!> \brief MO-native complex SEX: sigma_MO = prefactor * psi^T [ (psi rho psi^T) o W0 ] psi.
!>        Re/Im split through the real worker (complex-linear, P-C1 pattern).
!> \param bs_env ...
!> \param sigma_MO_cfm result, mo_active x mo_active (complex)
!> \param prefactor scaling applied to the final result
!> \param rho_MO_cfm input MO-active density-like matrix (complex)
!> \param ispin selects psi_sigma
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs_mo_complex(bs_env, sigma_MO_cfm, prefactor, rho_MO_cfm, ispin)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: sigma_MO_cfm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_MO_cfm
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs_mo_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_sigma

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%psi_built)

      CALL cp_fm_create(fm_rho, rho_MO_cfm%matrix_struct)
      CALL cp_fm_create(fm_sigma, rho_MO_cfm%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, rho_MO_cfm%matrix_struct)

      CALL cp_cfm_to_fm(msource=rho_MO_cfm, mtargetr=fm_rho)
      CALL sigma_ri_rs_worker(bs_env, fm_sigma, prefactor, fm_rho, &
                              bs_env%ri_rs%mat_psi_m_l(ispin))
      CALL cp_fm_to_cfm(msourcer=fm_sigma, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_MO_cfm, mtargeti=fm_rho)
      CALL sigma_ri_rs_worker(bs_env, fm_sigma, prefactor, fm_rho, &
                              bs_env%ri_rs%mat_psi_m_l(ispin))
      CALL cp_fm_to_cfm(msourcei=fm_sigma, mtarget=sigma_MO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), sigma_MO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_sigma)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs_mo_complex

! **************************************************************************************************
!> \brief MO-native real Hartree, single spin: V_MO = psi^T diag(v[n(rho)]) psi.
!> \param bs_env ...
!> \param rho_MO_fm input MO-active density (real)
!> \param V_MO_fm result, mo_active x mo_active (real)
!> \param ispin selects psi_sigma
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs_mo(bs_env, rho_MO_fm, V_MO_fm, ispin)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_MO_fm
      TYPE(cp_fm_type), INTENT(INOUT)                    :: V_MO_fm
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs_mo'

      INTEGER                                            :: handle, n_grid
      REAL(KIND=dp), ALLOCATABLE                         :: n_vec(:), v_vec(:)

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%psi_built)
      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, nfullrows_total=n_grid)
      ALLOCATE (n_vec(n_grid), v_vec(n_grid))
      n_vec(:) = 0.0_dp

      CALL hartree_density_ri_rs(bs_env, rho_MO_fm, bs_env%ri_rs%mat_psi_m_l(ispin), &
                                 n_vec, 1.0_dp)
      CALL hartree_potential_ri_rs(bs_env, n_vec, v_vec)
      CALL hartree_backtransform_ri_rs(bs_env, v_vec, V_MO_fm, &
                                       bs_env%ri_rs%mat_psi_m_l(ispin))

      DEALLOCATE (n_vec, v_vec)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs_mo

! **************************************************************************************************
!> \brief MO-native complex shared Hartree: one cross-spin potential v[sum_sigma g*n_sigma],
!>        per-spin back-transform V_MO_sigma = psi_sigma^T diag(v) psi_sigma. Re/Im split.
!> \param bs_env ...
!> \param rho_MO_cfms per-spin input MO-active densities (complex)
!> \param V_MO_cfms per-spin results (complex)
!> \param n_spin ...
!> \param spin_degeneracy density-accumulation factor g (2 closed shell, 1 open shell)
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs_mo_shared_complex(bs_env, rho_MO_cfms, V_MO_cfms, &
                                                      n_spin, spin_degeneracy)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), DIMENSION(:)                    :: rho_MO_cfms, V_MO_cfms
      INTEGER, INTENT(IN)                                :: n_spin
      REAL(KIND=dp), INTENT(IN)                          :: spin_degeneracy

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs_mo_shared_complex'

      INTEGER                                            :: handle, i_part, isp, n_grid
      REAL(KIND=dp), ALLOCATABLE                         :: n_vec(:), v_vec(:)
      TYPE(cp_cfm_type)                                  :: cfm_tmp
      TYPE(cp_fm_type)                                   :: fm_part, fm_v

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%psi_built)
      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL cp_fm_create(fm_part, rho_MO_cfms(1)%matrix_struct)
      CALL cp_fm_create(fm_v, rho_MO_cfms(1)%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, rho_MO_cfms(1)%matrix_struct)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, nfullrows_total=n_grid)
      ALLOCATE (n_vec(n_grid), v_vec(n_grid))

      DO i_part = 1, 2
         n_vec(:) = 0.0_dp
         DO isp = 1, n_spin
            IF (i_part == 1) THEN
               CALL cp_cfm_to_fm(msource=rho_MO_cfms(isp), mtargetr=fm_part)
            ELSE
               CALL cp_cfm_to_fm(msource=rho_MO_cfms(isp), mtargeti=fm_part)
            END IF
            CALL hartree_density_ri_rs(bs_env, fm_part, bs_env%ri_rs%mat_psi_m_l(isp), &
                                       n_vec, spin_degeneracy)
         END DO
         CALL hartree_potential_ri_rs(bs_env, n_vec, v_vec)
         DO isp = 1, n_spin
            CALL hartree_backtransform_ri_rs(bs_env, v_vec, fm_v, &
                                             bs_env%ri_rs%mat_psi_m_l(isp))
            IF (i_part == 1) THEN
               CALL cp_fm_to_cfm(msourcer=fm_v, mtarget=V_MO_cfms(isp))
            ELSE
               CALL cp_fm_to_cfm(msourcei=fm_v, mtarget=cfm_tmp)
               CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), V_MO_cfms(isp), &
                                         CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_tmp)
            END IF
         END DO
      END DO

      DEALLOCATE (n_vec, v_vec)
      CALL cp_fm_release(fm_part)
      CALL cp_fm_release(fm_v)
      CALL cp_cfm_release(cfm_tmp)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs_mo_shared_complex

! **************************************************************************************************
!> \brief Replicated row-wise dot: vec_out(l) = sum_n A(l,n)*B(l,n) for two dbcsr matrices on
!>        the same block structure/distribution. A-blocks without a B partner contribute zero.
!> \param matrix_A first factor (iterated)
!> \param matrix_B second factor (looked up locally; same distribution as A required)
!> \param vec_out full-length replicated result (zeroed on entry; sum-reduced on exit)
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE grid_row_dot(matrix_A, matrix_B, vec_out, para_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B
      REAL(KIND=dp), INTENT(INOUT)                       :: vec_out(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: col_blk, ic, ir, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_A, blk_B
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(matrix_A, nblkrows_total=nblkrows_total, row_blk_size=row_blk_size)
      ALLOCATE (row_offset(nblkrows_total + 1))
      row_offset(1) = 0
      DO ir = 1, nblkrows_total
         row_offset(ir + 1) = row_offset(ir) + row_blk_size(ir)
      END DO

      vec_out(:) = 0.0_dp

      CALL dbcsr_iterator_start(iter, matrix_A)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row_blk, col_blk, blk_A)
         CALL dbcsr_get_block_p(matrix_B, row_blk, col_blk, blk_B, found)
         IF (.NOT. found) CYCLE
         row_off = row_offset(row_blk)
         DO ic = 1, SIZE(blk_A, 2)
            DO ir = 1, SIZE(blk_A, 1)
               vec_out(row_off + ir) = vec_out(row_off + ir) + blk_A(ir, ic)*blk_B(ir, ic)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset)

      CALL para_env%sum(vec_out)

   END SUBROUTINE grid_row_dot

! **************************************************************************************************
!> \brief Scale each row of a dbcsr matrix by a replicated full-length vector:
!>        block(ir,ic) <- vec(global_row(ir)) * block(ir,ic). Value mutation only.
!> \param matrix ...
!> \param vec full-length replicated row-scaling vector
! **************************************************************************************************
   SUBROUTINE dbcsr_scale_rows_replicated(matrix, vec)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: vec(:)

      INTEGER                                            :: col_blk, ic, ir, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows_total, row_blk_size=row_blk_size)
      ALLOCATE (row_offset(nblkrows_total + 1))
      row_offset(1) = 0
      DO ir = 1, nblkrows_total
         row_offset(ir + 1) = row_offset(ir) + row_blk_size(ir)
      END DO

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row_blk, col_blk, blk)
         row_off = row_offset(row_blk)
         DO ic = 1, SIZE(blk, 2)
            DO ir = 1, SIZE(blk, 1)
               blk(ir, ic) = vec(row_off + ir)*blk(ir, ic)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset)

   END SUBROUTINE dbcsr_scale_rows_replicated

! **************************************************************************************************
!> \brief Replicated matvec: vec_out = matrix * vec_in (or matrix^T * vec_in if transposed) for a
!>        possibly rectangular distributed dbcsr matrix and replicated full-length vectors.
!>        Iterates over local blocks and reduces.
!> \param matrix distributed dbcsr matrix (may be rectangular)
!> \param vec_in full input vector, replicated on all ranks (column length, or row length if transposed)
!> \param vec_out full output vector, replicated on all ranks (zeroed on entry; sum-reduced on exit)
!> \param para_env ...
!> \param transposed if .TRUE. compute vec_out = matrix^T * vec_in
! **************************************************************************************************
   SUBROUTINE dbcsr_matvec_replicated(matrix, vec_in, vec_out, para_env, transposed)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: vec_in(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: vec_out(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: col_blk, col_off, ic, ir, &
                                                            nblkcols_total, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: col_offset(:), row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      LOGICAL                                            :: my_trans
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      my_trans = .FALSE.
      IF (PRESENT(transposed)) my_trans = transposed

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)

      ALLOCATE (row_offset(nblkrows_total + 1), col_offset(nblkcols_total + 1))
      row_offset(1) = 0
      DO ir = 1, nblkrows_total
         row_offset(ir + 1) = row_offset(ir) + row_blk_size(ir)
      END DO
      col_offset(1) = 0
      DO ic = 1, nblkcols_total
         col_offset(ic + 1) = col_offset(ic) + col_blk_size(ic)
      END DO

      vec_out(:) = 0.0_dp

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row_blk, col_blk, block)
         row_off = row_offset(row_blk)
         col_off = col_offset(col_blk)
         IF (my_trans) THEN
            DO ic = 1, SIZE(block, 2)
               DO ir = 1, SIZE(block, 1)
                  vec_out(col_off + ic) = vec_out(col_off + ic) + &
                                          block(ir, ic)*vec_in(row_off + ir)
               END DO
            END DO
         ELSE
            DO ic = 1, SIZE(block, 2)
               DO ir = 1, SIZE(block, 1)
                  vec_out(row_off + ir) = vec_out(row_off + ir) + &
                                          block(ir, ic)*vec_in(col_off + ic)
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset, col_offset)

      CALL para_env%sum(vec_out)

   END SUBROUTINE dbcsr_matvec_replicated

END MODULE rt_bse_ri_rs
