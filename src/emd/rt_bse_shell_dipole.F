!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distance-resolved linRTBSE diagnostics on the AO density, sharing one binning geometry.
!>        (i) the printed MOMENTS trace split by the inter-atomic distance of the AO pair that
!>        carries it,
!>          µ^ν_b(t) = -g sum_µλ in b  dT^ν_µλ ρ^AO_µλ,   b = bin(|R_A(µ) - R_B(λ)|)
!>        so that sum_b µ^ν_b(t) = µ^ν(t) exactly; (ii) the running Fourier transform of the AO
!>        density at selected frequencies Ω, binned by the same distance at end of run, giving the
!>        distance profile of the response at Ω rather than of its dipole alone. Own module because
!>        the binning geometry comes from bs_env%ri_rs while the density comes from rtbse_env.
!>        [INTERNAL] dev_notes/rt_bse_freq_resolved_distance_handoff.md
!> \author Max Graml
! **************************************************************************************************
MODULE rt_bse_shell_dipole

   USE block_decay_diag,                ONLY: build_block_decay_maps
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_printkey_is_on
   USE gw_non_periodic_ri_rs,           ONLY: get_basis_offsets
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              femtoseconds
   USE qs_environment_types,            ONLY: get_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE rt_bse_ri_rs,                    ONLY: compute_sigma_ri_rs_complex,&
                                              sex_energy_by_distance
   USE rt_bse_types,                    ONLY: rtbse_env_type
   USE rt_propagation_ft,               ONLY: multi_fft
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_shell_dipole'

   PUBLIC :: accumulate_shell_dipole, &
             accumulate_shell_rho_ft, &
             build_shell_dipole_maps, &
             print_shell_dipole, &
             print_shell_dipole_ft, &
             print_shell_rho_ft, &
             print_shell_sex_ft

CONTAINS

! **************************************************************************************************
!> \brief Build what the distance-resolved diagnostics need: the distance bins (shared with the
!>        block-decay diagnostic, so the two are comparable bin-for-bin), the AO -> atom map, the
!>        AO-sized scratch, the dipole trace and, when SHELL_RHO_OMEGA is set, the ρ̃(Ω)
!>        accumulators. Must run after create_rtbse_env (which reads the keywords and
!>        resolves the RI-RS flag) and before initialize_maximum_timestep, whose
!>        reallocate_ft_traces resizes the trace if ENFORCE_MAX_DT changes the step count, and
!>        before initialize_moments, which retains the AO operators only when the flag is set.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
! **************************************************************************************************
   SUBROUTINE build_shell_dipole_maps(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_shell_dipole_maps'

      INTEGER                                            :: handle, iao, iatom, isp, iw, n_omega, &
                                                            natom, total_sgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      ! do_shell_dipole comes from &RTBSE SHELL_DIPOLE, not from the print keys: cp_printkey_is_on
      ! is global_print_level >= key_print_level, so a print key cannot express "off by default"
      IF (.NOT. rtbse_env%do_shell_dipole) RETURN

      CALL timeset(routineN, handle)

      ! The bins are the block-decay bins, and those live on the RI-RS grid caches
      IF (.NOT. rtbse_env%rirs_kernel) THEN
         CALL cp_abort(__LOCATION__, &
                       "SHELL_DIPOLE needs the real-space kernel: set RTBSE%KERNEL_RI to RS "// &
                       "and activate the &GW_RI_RS section, or unset RTBSE%SHELL_DIPOLE.")
      END IF
      IF (.NOT. ALLOCATED(rtbse_env%bs_env%ri_rs%atom_centers) .OR. &
          .NOT. ALLOCATED(rtbse_env%bs_env%ri_rs%chunk_centroids)) THEN
         CALL cp_abort(__LOCATION__, &
                       "SHELL_DIPOLE: the RI-RS grid caches are absent, so the distance bins "// &
                       "cannot be built. Activate the &GW_RI_RS section.")
      END IF
      CALL build_block_decay_maps(rtbse_env%bs_env)
      rtbse_env%shell_nbins = rtbse_env%bs_env%ri_rs%decay_nbins

      ! AO index -> atom index, so an AO pair inherits the distance of its parent atoms
      CALL get_qs_env(rtbse_env%qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
      natom = SIZE(particle_set)
      ALLOCATE (first_sgf(natom + 1))
      CALL get_basis_offsets(particle_set, qs_kind_set, first_sgf, total_sgf)
      CPASSERT(total_sgf == rtbse_env%n_ao)
      ALLOCATE (rtbse_env%shell_ao_atom(total_sgf))
      DO iatom = 1, natom
         DO iao = first_sgf(iatom), first_sgf(iatom + 1) - 1
            rtbse_env%shell_ao_atom(iao) = iatom
         END DO
      END DO
      DEALLOCATE (first_sgf)

      ! Same struct as moments_ao, so local element (ii,jj) is the same global pair in both
      CALL cp_cfm_create(rtbse_env%rho_ao_shell, rtbse_env%bs_env%fm_ks_Gamma(1)%matrix_struct)

      ! Primary allocation, mirroring moments_trace in create_rtbse_env; reallocate_ft_traces only
      ! resizes it, and only when ENFORCE_MAX_DT changes the step count
      ALLOCATE (rtbse_env%shell_dipole_trace(rtbse_env%n_spin, 3, rtbse_env%shell_nbins, &
                                             rtbse_env%sim_nsteps + 1), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      ! ρ̃(Ω) accumulators: one AO matrix per requested frequency and spin, plus the sample that
      ! multi_fft's subtract_initial would remove. Distributed, so the footprint below is the total.
      IF (rtbse_env%do_shell_rho_ft) THEN
         n_omega = SIZE(rtbse_env%shell_rho_omega)
         ALLOCATE (rtbse_env%shell_rho_ref(rtbse_env%n_spin))
         ALLOCATE (rtbse_env%shell_rho_ft(n_omega, rtbse_env%n_spin))
         DO isp = 1, rtbse_env%n_spin
            CALL cp_cfm_create(rtbse_env%shell_rho_ref(isp), rtbse_env%bs_env%fm_ks_Gamma(1)%matrix_struct)
            CALL cp_cfm_set_all(rtbse_env%shell_rho_ref(isp), CMPLX(0.0_dp, 0.0_dp, KIND=dp))
            DO iw = 1, n_omega
               CALL cp_cfm_create(rtbse_env%shell_rho_ft(iw, isp), rtbse_env%bs_env%fm_ks_Gamma(1)%matrix_struct)
               CALL cp_cfm_set_all(rtbse_env%shell_rho_ft(iw, isp), CMPLX(0.0_dp, 0.0_dp, KIND=dp))
            END DO
         END DO
         IF (rtbse_env%unit_nr > 0) THEN
            WRITE (rtbse_env%unit_nr, '(T2,A,T61,I20)') "SHELL_RHO_FT| Number of frequencies", n_omega
            WRITE (rtbse_env%unit_nr, '(T2,A,T61,F20.1)') "SHELL_RHO_FT| Matrix storage, all ranks [MB]", &
               REAL(n_omega + 1, dp)*REAL(rtbse_env%n_spin, dp)* &
               (REAL(rtbse_env%n_ao, dp)**2)*16.0_dp/1.0E6_dp
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_shell_dipole_maps

! **************************************************************************************************
!> \brief Bin one spin channel's dipole contribution by AO-pair distance:
!>          µ^ν_b = -g sum_µλ in b  dT^ν_µλ ρ^AO_µλ,   b = bin(|R_A(µ) - R_B(λ)|)
!>        dT^ν = (d^ν)^T is stored transposed at setup, so this is the printed MOMENTS integrand
!>        -g sum_µλ d^ν_λµ ρ^AO_µλ element for element, with no symmetry assumption on d^ν.
!>        Rank-local over the ScaLAPACK block-cyclic tiles, then summed - the loop touches only
!>        local_data, never a collective element accessor.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param rho_ao Lab-frame ρ^AO for this spin, from transform_mo_to_ao_contravariant_cfm
!> \param isp Spin channel σ
! **************************************************************************************************
   SUBROUTINE accumulate_shell_dipole(rtbse_env, rho_ao, isp)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho_ao
      INTEGER, INTENT(IN)                                :: isp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'accumulate_shell_dipole'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: acc
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: rho_local
      INTEGER                                            :: handle, ibin, ii, jj, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: r_ij
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: d_x, d_y, d_z
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1
      ALLOCATE (acc(3, rtbse_env%shell_nbins), source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      CALL cp_cfm_get_info(rho_ao, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices, &
                           local_data=rho_local, para_env=para_env)
      CALL cp_fm_get_info(rtbse_env%moments_ao(1), local_data=d_x)
      CALL cp_fm_get_info(rtbse_env%moments_ao(2), local_data=d_y)
      CALL cp_fm_get_info(rtbse_env%moments_ao(3), local_data=d_z)

      ! bin index depends on the pair only, so it is hoisted out of the Cartesian triple
      DO jj = 1, ncol_local
         DO ii = 1, nrow_local
            r_ij = NORM2(rtbse_env%bs_env%ri_rs%atom_centers(:, rtbse_env%shell_ao_atom(row_indices(ii))) &
                         - rtbse_env%bs_env%ri_rs%atom_centers(:, rtbse_env%shell_ao_atom(col_indices(jj))))
            ibin = MIN(INT(r_ij/rtbse_env%bs_env%ri_rs%decay_bin_width) + 1, rtbse_env%shell_nbins)
            acc(1, ibin) = acc(1, ibin) + d_x(ii, jj)*rho_local(ii, jj)
            acc(2, ibin) = acc(2, ibin) + d_y(ii, jj)*rho_local(ii, jj)
            acc(3, ibin) = acc(3, ibin) + d_z(ii, jj)*rho_local(ii, jj)
         END DO
      END DO

      ! collective over the matrix's own communicator, outside every rank-conditional branch
      CALL para_env%sum(acc)

      DO ibin = 1, rtbse_env%shell_nbins
         DO k = 1, 3
            rtbse_env%shell_dipole_trace(isp, k, ibin, n) = -acc(k, ibin)*rtbse_env%spin_degeneracy
         END DO
      END DO

      DEALLOCATE (acc)
      CALL timestop(handle)

   END SUBROUTINE accumulate_shell_dipole

! **************************************************************************************************
!> \brief Emit one block of rows per step: one row per populated distance bin, per spin. The first
!>        recorded step is the pre-kick static profile - it is what multi_fft's subtract_initial
!>        removes, and it is flagged in the file rather than written to a separate channel.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param t_phys Physical time associated with the recorded density
!> \param n Time index into shell_dipole_trace, shared with moments_trace
! **************************************************************************************************
   SUBROUTINE print_shell_dipole(rtbse_env, t_phys, n)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      REAL(KIND=dp), INTENT(IN)                          :: t_phys
      INTEGER, INTENT(IN)                                :: n

      CHARACTER(LEN=11), DIMENSION(2), PARAMETER :: &
         file_extensions = ['_SPIN_A.dat', '_SPIN_B.dat']

      INTEGER                                            :: ibin, isp, print_unit
      REAL(KIND=dp)                                      :: r_ang
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      IF (.NOT. cp_printkey_is_on(logger%iter_info, rtbse_env%shell_dipole_section)) RETURN

      DO isp = 1, rtbse_env%n_spin
         print_unit = cp_print_key_unit_nr(logger, rtbse_env%shell_dipole_section, &
                                           extension=file_extensions(isp))
         IF (print_unit > 0) THEN
            IF (.NOT. rtbse_env%shell_hdr_written) THEN
               ! full mantissa: the sum rule against MOMENTS is gated at round-off, and a
               ! 7-significant-digit column would cap the achievable agreement at ~1e-7
               WRITE (print_unit, '(A)') "#      Time [fs]  bin   R [Ang]"// &
                  "                    x,re                    x,im                    y,re"// &
                  "                    y,im                    z,re                    z,im"
               IF (rtbse_env%restart_extracted) THEN
                  WRITE (print_unit, '(A)') "# RESTART: rows before this run's first step are "// &
                     "absent from the trace, so no spectrum is emitted - see SHELL_DIPOLE_FT"
               ELSE
                  WRITE (print_unit, '(A)') "# the block below is the pre-kick static reference profile"
               END IF
            END IF
            DO ibin = 1, rtbse_env%shell_nbins
               IF (rtbse_env%bs_env%ri_rs%decay_npairs(ibin, 1) == 0) CYCLE
               r_ang = (REAL(ibin, dp) - 0.5_dp)*rtbse_env%bs_env%ri_rs%decay_bin_width*angstrom
               WRITE (print_unit, '(F16.6,I5,F10.2,6ES24.15)') &
                  t_phys*femtoseconds, ibin, r_ang, &
                  REAL(rtbse_env%shell_dipole_trace(isp, 1, ibin, n), KIND=dp), &
                  AIMAG(rtbse_env%shell_dipole_trace(isp, 1, ibin, n)), &
                  REAL(rtbse_env%shell_dipole_trace(isp, 2, ibin, n), KIND=dp), &
                  AIMAG(rtbse_env%shell_dipole_trace(isp, 2, ibin, n)), &
                  REAL(rtbse_env%shell_dipole_trace(isp, 3, ibin, n), KIND=dp), &
                  AIMAG(rtbse_env%shell_dipole_trace(isp, 3, ibin, n))
            END DO
         END IF
         CALL cp_print_key_finished_output(print_unit, logger, rtbse_env%shell_dipole_section)
      END DO
      ! after the spin loop, so both per-spin files carry the header
      rtbse_env%shell_hdr_written = .TRUE.

   END SUBROUTINE print_shell_dipole

! **************************************************************************************************
!> \brief End of run: transform every per-bin dipole trace,
!>          µ^ν_b(ω) = dt sum_t [µ^ν_b(t) - µ^ν_b(t_0)] exp(-γ(t-t_0)) exp(-i ω (t-t_0))
!>        through the SAME multi_fft call the MOMENTS_FT path uses, with the same damping γ,
!>        reference time t_0 and initial subtraction. sum_b µ^ν_b(ω) = µ^ν(ω) then holds by
!>        linearity of the transform, so closure against MOMENTS_FT is exact rather than a
!>        re-implemented window agreeing to plotting accuracy.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
! **************************************************************************************************
   SUBROUTINE print_shell_dipole_ft(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_shell_dipole_ft'
      CHARACTER(LEN=11), DIMENSION(2), PARAMETER :: &
         file_extensions = ['_SPIN_A.dat', '_SPIN_B.dat']

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: results, series
      INTEGER                                            :: handle, ibin, iser, isp, iw, k, n, nser, &
                                                            print_unit
      REAL(KIND=dp)                                      :: r_ang
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omegas
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (.NOT. rtbse_env%do_shell_dipole) RETURN
      logger => cp_get_default_logger()
      IF (.NOT. cp_printkey_is_on(logger%iter_info, rtbse_env%shell_dipole_ft_section)) RETURN

      ! read_restart_trace restores moments_trace but has no shell counterpart, so on a restart the
      ! slots before this run's first step are still zero. Transforming that gives a spectrum with a
      ! zero prefix and subtract_initial removes 0 instead of the reference - refuse rather than
      ! emit a file whose header claims closure against MOMENTS_FT.
      IF (rtbse_env%restart_extracted) THEN
         CALL cp_warn(__LOCATION__, &
                      "SHELL_DIPOLE_FT skipped on a restart: the shell trace is not restored, so "// &
                      "its spectrum would not close on MOMENTS_FT. The time-domain rows are still written.")
         RETURN
      END IF

      CALL timeset(routineN, handle)

      n = SIZE(rtbse_env%time_trace)
      nser = 3*rtbse_env%shell_nbins
      ALLOCATE (series(nser, n))
      ALLOCATE (results(nser, n))
      ALLOCATE (omegas(n))

      DO isp = 1, rtbse_env%n_spin
         DO ibin = 1, rtbse_env%shell_nbins
            DO k = 1, 3
               series(3*(ibin - 1) + k, :) = rtbse_env%shell_dipole_trace(isp, k, ibin, :)
            END DO
         END DO
         CALL multi_fft(rtbse_env%time_trace, series, results, omegas, &
                        damping_opt=rtbse_env%dft_control%rtp_control%ft_damping, &
                        t0_opt=rtbse_env%dft_control%rtp_control%ft_t0, &
                        subtract_initial_opt=.TRUE.)
         print_unit = cp_print_key_unit_nr(logger, rtbse_env%shell_dipole_ft_section, &
                                           extension=file_extensions(isp), &
                                           file_form="FORMATTED", file_position="REWIND")
         IF (print_unit > 0) THEN
            ! omega in a.u., matching the MOMENTS_FT / POLARIZABILITY .dat convention - the
            ! closure check joins the two files on this column
            WRITE (print_unit, '(A)') "#  omega [at.u.]  bin   R [Ang]"// &
               "                    x,re                    x,im                    y,re"// &
               "                    y,im                    z,re                    z,im"
            ! multi_fft truncates to the samples at or after t0, so SIZE(results,2) - not n -
            ! is the number of frequencies it actually produced
            DO iw = 1, SIZE(results, 2)
               DO ibin = 1, rtbse_env%shell_nbins
                  IF (rtbse_env%bs_env%ri_rs%decay_npairs(ibin, 1) == 0) CYCLE
                  r_ang = (REAL(ibin, dp) - 0.5_dp)*rtbse_env%bs_env%ri_rs%decay_bin_width*angstrom
                  iser = 3*(ibin - 1)
                  WRITE (print_unit, '(F16.6,I5,F10.2,6ES24.15)') &
                     omegas(iw), ibin, r_ang, &
                     (REAL(results(iser + k, iw), KIND=dp), AIMAG(results(iser + k, iw)), k=1, 3)
               END DO
            END DO
         END IF
         CALL cp_print_key_finished_output(print_unit, logger, rtbse_env%shell_dipole_ft_section)
      END DO

      DEALLOCATE (series, results, omegas)
      CALL timestop(handle)

   END SUBROUTINE print_shell_dipole_ft

! **************************************************************************************************
!> \brief Add one time sample to the running Fourier transform of the AO density,
!>          ρ̃_µλ(Ω) += dt [ρ^AO_µλ(t) - ρ^AO_µλ(t_1)] exp(-γ (t - t_ref)) exp(+i Ω (t - t_ref))
!>        which is multi_fft's kernel written as a sum over steps: same dt scaling, same damping γ,
!>        same subtraction of the first recorded sample. The sign is multi_fft's measured one -
!>        exp(+iΩt), not the exp(-iΩt) the fsign=-1 plan argument suggests. Ω is taken literally
!>        from SHELL_RHO_OMEGA and means 2π/period, whereas multi_fft LABELS its bin m with
!>        2πm/((N-1)dt) while transforming it at 2πm/(N dt) - the two axes differ by N/(N-1).
!>        The accumulation is elementwise on the rank-local tiles - no communication per step.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param rho_ao Lab-frame ρ^AO for this spin, from transform_mo_to_ao_contravariant_cfm
!> \param isp Spin channel σ
!> \param t_phys Physical time of this sample
!> \param n Time index into the traces, shared with moments_trace
! **************************************************************************************************
   SUBROUTINE accumulate_shell_rho_ft(rtbse_env, rho_ao, isp, t_phys, n)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho_ao
      INTEGER, INTENT(IN)                                :: isp
      REAL(KIND=dp), INTENT(IN)                          :: t_phys
      INTEGER, INTENT(IN)                                :: n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'accumulate_shell_rho_ft'

      COMPLEX(KIND=dp)                                   :: zfac
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: ft_local, ref_local, rho_local
      INTEGER                                            :: handle, iw
      REAL(KIND=dp)                                      :: damping_opt, dt_rel, weight

      IF (.NOT. rtbse_env%do_shell_rho_ft) RETURN
      ! the reference sample and the window origin are both taken from this run's own trace, which
      ! a restart does not carry - print_shell_rho_ft refuses to emit for the same reason
      IF (rtbse_env%restart_extracted) RETURN

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(rho_ao, local_data=rho_local)
      CALL cp_cfm_get_info(rtbse_env%shell_rho_ref(isp), local_data=ref_local)

      ! multi_fft subtracts value_series(:,1), the first recorded sample, whatever the window start
      IF (n == 1) ref_local(:, :) = rho_local(:, :)

      IF (t_phys < rtbse_env%dft_control%rtp_control%ft_t0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (rtbse_env%shell_ft_tref < 0.0_dp) THEN
         rtbse_env%shell_ft_tref = t_phys
         ! γ = 4 / (window length) unless FT%DAMPING gives a decay time; mirrors multi_fft exactly,
         ! and the window end is known here because ENFORCE_MAX_DT has already fixed sim_nsteps
         damping_opt = rtbse_env%dft_control%rtp_control%ft_damping
         IF (damping_opt > 0.0_dp) THEN
            rtbse_env%shell_ft_gamma = 1.0_dp/damping_opt
         ELSE IF (damping_opt == 0.0_dp) THEN
            rtbse_env%shell_ft_gamma = 0.0_dp
         ELSE
            rtbse_env%shell_ft_gamma = 4.0_dp/ &
                                       (REAL(rtbse_env%sim_nsteps, dp)*rtbse_env%sim_dt - t_phys)
         END IF
      END IF

      IF (isp == 1) rtbse_env%shell_ft_nsamples = rtbse_env%shell_ft_nsamples + 1

      dt_rel = t_phys - rtbse_env%shell_ft_tref
      weight = rtbse_env%sim_dt*EXP(-rtbse_env%shell_ft_gamma*dt_rel)

      DO iw = 1, SIZE(rtbse_env%shell_rho_omega)
         zfac = weight*EXP(CMPLX(0.0_dp, rtbse_env%shell_rho_omega(iw)*dt_rel, KIND=dp))
         CALL cp_cfm_get_info(rtbse_env%shell_rho_ft(iw, isp), local_data=ft_local)
         ft_local(:, :) = ft_local(:, :) + zfac*(rho_local(:, :) - ref_local(:, :))
      END DO

      CALL timestop(handle)

   END SUBROUTINE accumulate_shell_rho_ft

! **************************************************************************************************
!> \brief End of run: bin ρ̃(Ω) by the same AO-pair distance the dipole uses and emit, per frequency
!>        and populated bin,
!>          v_b(Ω) = sqrt( sum_µλ in b |ρ̃_µλ(Ω)|² )   and   µ^ν_b(Ω) = -g sum_µλ in b dT^ν_µλ ρ̃_µλ(Ω)
!>        v_b is the distance profile of the response itself; µ^ν_b(Ω) is the same object the FT of
!>        the shell dipole produces bin for bin, so the two files cross-check the two orders of
!>        (bin, transform). Join them by FFT bin index, m = round(Ω N dt / 2π), never by the
!>        printed ω - the two files use different ω axes, see accumulate_shell_rho_ft.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
! **************************************************************************************************
   SUBROUTINE print_shell_rho_ft(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_shell_rho_ft'
      CHARACTER(LEN=11), DIMENSION(2), PARAMETER :: &
         file_extensions = ['_SPIN_A.dat', '_SPIN_B.dat']

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: dipole
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: ft_local
      INTEGER                                            :: handle, ibin, ii, isp, iw, jj, k, &
                                                            ncol_local, nrow_local, print_unit
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: r_ang, r_ij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: profile
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: d_x, d_y, d_z
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env

      IF (.NOT. rtbse_env%do_shell_rho_ft) RETURN
      logger => cp_get_default_logger()
      IF (.NOT. cp_printkey_is_on(logger%iter_info, rtbse_env%shell_rho_ft_section)) RETURN

      IF (rtbse_env%restart_extracted) THEN
         CALL cp_warn(__LOCATION__, &
                      "SHELL_RHO_FT skipped on a restart: the steps before the continuation are "// &
                      "absent from the accumulator, so its spectrum would not match SHELL_DIPOLE_FT.")
         RETURN
      END IF

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(rtbse_env%rho_ao_shell, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices, para_env=para_env)
      CALL cp_fm_get_info(rtbse_env%moments_ao(1), local_data=d_x)
      CALL cp_fm_get_info(rtbse_env%moments_ao(2), local_data=d_y)
      CALL cp_fm_get_info(rtbse_env%moments_ao(3), local_data=d_z)

      ALLOCATE (profile(rtbse_env%shell_nbins, SIZE(rtbse_env%shell_rho_omega)), source=0.0_dp)
      ALLOCATE (dipole(3, rtbse_env%shell_nbins, SIZE(rtbse_env%shell_rho_omega)), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      DO isp = 1, rtbse_env%n_spin
         profile(:, :) = 0.0_dp
         dipole(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO iw = 1, SIZE(rtbse_env%shell_rho_omega)
            CALL cp_cfm_get_info(rtbse_env%shell_rho_ft(iw, isp), local_data=ft_local)
            DO jj = 1, ncol_local
               DO ii = 1, nrow_local
                  r_ij = NORM2(rtbse_env%bs_env%ri_rs%atom_centers(:, rtbse_env%shell_ao_atom(row_indices(ii))) &
                               - rtbse_env%bs_env%ri_rs%atom_centers(:, rtbse_env%shell_ao_atom(col_indices(jj))))
                  ibin = MIN(INT(r_ij/rtbse_env%bs_env%ri_rs%decay_bin_width) + 1, rtbse_env%shell_nbins)
                  profile(ibin, iw) = profile(ibin, iw) + ABS(ft_local(ii, jj))**2
                  dipole(1, ibin, iw) = dipole(1, ibin, iw) + d_x(ii, jj)*ft_local(ii, jj)
                  dipole(2, ibin, iw) = dipole(2, ibin, iw) + d_y(ii, jj)*ft_local(ii, jj)
                  dipole(3, ibin, iw) = dipole(3, ibin, iw) + d_z(ii, jj)*ft_local(ii, jj)
               END DO
            END DO
         END DO
         ! collective, outside every rank-conditional branch
         CALL para_env%sum(profile)
         CALL para_env%sum(dipole)
         dipole(:, :, :) = -dipole(:, :, :)*rtbse_env%spin_degeneracy

         print_unit = cp_print_key_unit_nr(logger, rtbse_env%shell_rho_ft_section, &
                                           extension=file_extensions(isp), &
                                           file_form="FORMATTED", file_position="REWIND")
         IF (print_unit > 0) THEN
            ! the window is recorded in full precision so the file can be checked against the
            ! SHELL_DIPOLE time trace without reconstructing dt from a rounded time column
            WRITE (print_unit, '(A,3ES24.15,A,I0)') "# window t_ref, dt, gamma [at.u.]:", &
               rtbse_env%shell_ft_tref, rtbse_env%sim_dt, rtbse_env%shell_ft_gamma, &
               "   nsamples: ", rtbse_env%shell_ft_nsamples
            WRITE (print_unit, '(A)') "# omega in a.u., kernel exp(+i w t); |rho| carries no spin"// &
               " degeneracy factor, the dipole columns do"
            WRITE (print_unit, '(A)') "#            omega [at.u.]  bin   R [Ang]"// &
               "                  |rho|                    x,re                    x,im"// &
               "                    y,re                    y,im                    z,re"// &
               "                    z,im"
            DO iw = 1, SIZE(rtbse_env%shell_rho_omega)
               DO ibin = 1, rtbse_env%shell_nbins
                  IF (rtbse_env%bs_env%ri_rs%decay_npairs(ibin, 1) == 0) CYCLE
                  r_ang = (REAL(ibin, dp) - 0.5_dp)*rtbse_env%bs_env%ri_rs%decay_bin_width*angstrom
                  WRITE (print_unit, '(ES24.15,I5,F10.2,7ES24.15)') &
                     rtbse_env%shell_rho_omega(iw), ibin, r_ang, SQRT(profile(ibin, iw)), &
                     (REAL(dipole(k, ibin, iw), KIND=dp), AIMAG(dipole(k, ibin, iw)), k=1, 3)
               END DO
            END DO
         END IF
         CALL cp_print_key_finished_output(print_unit, logger, rtbse_env%shell_rho_ft_section)
      END DO

      DEALLOCATE (profile, dipole)
      CALL timestop(handle)

   END SUBROUTINE print_shell_rho_ft

! **************************************************************************************************
!> \brief End of run: for each stored frequency, split the SEX energy of ρ̃(Ω) over the W^0 block
!>        distance,
!>          E_b(Ω) = sum_ll' in b  W^0_ll' |ρ̃^grid_ll'(Ω)|²,   sum_b E_b(Ω) = Tr[ρ̃^† Σ^SEX[ρ̃]]
!>        and emit it with its running cumulative sum, so "truncating W^0 at R keeps this fraction
!>        of the response at Ω" is read straight off the file. The closure value on each block
!>        header is computed independently, through the production kernel
!>        compute_sigma_ri_rs_complex - it is the gate, not a restatement.
!>        ⚠ These bins are GRID-chunk distances; SHELL_RHO_FT's are AO parent-atom distances.
!>        Same widths, different populations - do not join the two files on the bin index.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
! **************************************************************************************************
   SUBROUTINE print_shell_sex_ft(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_shell_sex_ft'
      CHARACTER(LEN=11), DIMENSION(2), PARAMETER :: &
         file_extensions = ['_SPIN_A.dat', '_SPIN_B.dat']

      COMPLEX(KIND=dp)                                   :: trace_sex
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: rho_local, sigma_local
      INTEGER                                            :: handle, ibin, ic, isp, iw, jc, nbins, &
                                                            print_unit
      REAL(KIND=dp)                                      :: cum_e, r_ang, r_bins, r_chunk_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_bin, s_bin
      TYPE(cp_cfm_type)                                  :: sigma_cfm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env

      IF (.NOT. rtbse_env%do_shell_rho_ft) RETURN
      logger => cp_get_default_logger()
      IF (.NOT. cp_printkey_is_on(logger%iter_info, rtbse_env%shell_sex_ft_section)) RETURN

      IF (rtbse_env%restart_extracted) THEN
         CALL cp_warn(__LOCATION__, &
                      "SHELL_SEX_FT skipped on a restart: it is built from the same accumulator "// &
                      "as SHELL_RHO_FT, which the restart does not carry.")
         RETURN
      END IF
      IF (.NOT. rtbse_env%bs_env%ri_rs%W0_grid_built) THEN
         CALL cp_warn(__LOCATION__, &
                      "SHELL_SEX_FT skipped: the RI-RS W^0 grid was never built, so there is no "// &
                      "SEX kernel to split. Unset RTBSE%DEBUG_DISABLE_SEX.")
         RETURN
      END IF

      CALL timeset(routineN, handle)

      nbins = rtbse_env%bs_env%ri_rs%decay_nbins
      ALLOCATE (e_bin(nbins), s_bin(nbins))

      ! The bin count comes from the ATOM bounding box, but these are GRID distances, which reach
      ! further - so the last bin is a clamp holding every pair beyond it. That matters here more
      ! than for the max-norm slots: the far tail is the quantity this file exists to measure.
      r_chunk_max = 0.0_dp
      DO ic = 1, SIZE(rtbse_env%bs_env%ri_rs%chunk_centroids, 2)
         DO jc = ic + 1, SIZE(rtbse_env%bs_env%ri_rs%chunk_centroids, 2)
            r_chunk_max = MAX(r_chunk_max, NORM2(rtbse_env%bs_env%ri_rs%chunk_centroids(:, ic) - &
                                                 rtbse_env%bs_env%ri_rs%chunk_centroids(:, jc)))
         END DO
      END DO
      r_bins = REAL(nbins - 1, dp)*rtbse_env%bs_env%ri_rs%decay_bin_width
      CALL cp_cfm_create(sigma_cfm, rtbse_env%bs_env%fm_ks_Gamma(1)%matrix_struct)

      DO isp = 1, rtbse_env%n_spin
         print_unit = cp_print_key_unit_nr(logger, rtbse_env%shell_sex_ft_section, &
                                           extension=file_extensions(isp), &
                                           file_form="FORMATTED", file_position="REWIND")
         IF (print_unit > 0) THEN
            WRITE (print_unit, '(A)') "# SEX energy of the response at omega, split by the W^0 "// &
               "block distance |r_l - r_l'| on the grid chunk centroids."
            WRITE (print_unit, '(A)') "# Bins are GRID distances - NOT the AO parent-atom bins of "// &
               "SHELL_RHO_FT. cum_E_b is the running sum over bins."
            WRITE (print_unit, '(A,I0,A,F8.2,A,F8.2,A)') &
               "# last bin (", nbins, ") is a CLAMP: it holds every pair beyond ", r_bins*angstrom, &
               " Ang, and grid chunks reach ", r_chunk_max*angstrom, " Ang apart"
            WRITE (print_unit, '(A)') "#            omega [at.u.]  bin   R [Ang]"// &
               "                    E_b                cum_E_b                    S_b"
         END IF
         DO iw = 1, SIZE(rtbse_env%shell_rho_omega)
            CALL sex_energy_by_distance(rtbse_env%bs_env, rtbse_env%shell_rho_ft(iw, isp), &
                                        e_bin, s_bin)
            ! independent closure: the same energy through the production kernel, at prefactor 1
            CALL compute_sigma_ri_rs_complex(rtbse_env%bs_env, sigma_cfm, 1.0_dp, &
                                             rtbse_env%shell_rho_ft(iw, isp))
            CALL cp_cfm_get_info(rtbse_env%shell_rho_ft(iw, isp), local_data=rho_local, &
                                 para_env=para_env)
            CALL cp_cfm_get_info(sigma_cfm, local_data=sigma_local)
            trace_sex = SUM(CONJG(rho_local(:, :))*sigma_local(:, :))
            CALL para_env%sum(trace_sex)

            IF (print_unit > 0) THEN
               WRITE (print_unit, '(A,ES24.15,A,ES24.15,A,2ES24.15)') &
                  "# omega ", rtbse_env%shell_rho_omega(iw), "   sum_b E_b ", SUM(e_bin), &
                  "   Tr[rho~^dag Sigma] ", REAL(trace_sex, KIND=dp), AIMAG(trace_sex)
               cum_e = 0.0_dp
               DO ibin = 1, nbins
                  cum_e = cum_e + e_bin(ibin)
                  IF (rtbse_env%bs_env%ri_rs%decay_npairs(ibin, 2) == 0) CYCLE
                  r_ang = (REAL(ibin, dp) - 0.5_dp)*rtbse_env%bs_env%ri_rs%decay_bin_width*angstrom
                  WRITE (print_unit, '(ES24.15,I5,F10.2,3ES24.15)') &
                     rtbse_env%shell_rho_omega(iw), ibin, r_ang, e_bin(ibin), cum_e, s_bin(ibin)
               END DO
            END IF
         END DO
         CALL cp_print_key_finished_output(print_unit, logger, rtbse_env%shell_sex_ft_section)
      END DO

      CALL cp_cfm_release(sigma_cfm)
      DEALLOCATE (e_bin, s_bin)
      CALL timestop(handle)

   END SUBROUTINE print_shell_sex_ft

END MODULE rt_bse_shell_dipole
