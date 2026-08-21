!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shell-resolved dipole for the linearized RT-BSE: the printed MOMENTS trace split by the
!>        inter-atomic distance of the AO pair that carries it,
!>          µ^ν_b(t) = -g sum_µλ in b  dT^ν_µλ ρ^AO_µλ,   b = bin(|R_A(µ) - R_B(λ)|)
!>        so that sum_b µ^ν_b(t) = µ^ν(t) exactly. Own module because the binning geometry comes
!>        from bs_env%ri_rs while the density comes from rtbse_env.
!>        [INTERNAL] dev_notes/rt_bse_freq_resolved_distance_handoff.md
!> \author Max Graml
! **************************************************************************************************
MODULE rt_bse_shell_dipole

   USE block_decay_diag,                ONLY: build_block_decay_maps
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
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
   USE rt_bse_types,                    ONLY: rtbse_env_type
   USE rt_propagation_ft,               ONLY: multi_fft
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_shell_dipole'

   PUBLIC :: accumulate_shell_dipole, &
             build_shell_dipole_maps, &
             print_shell_dipole, &
             print_shell_dipole_ft

CONTAINS

! **************************************************************************************************
!> \brief Build what the shell dipole needs: the distance bins (shared with the block-decay
!>        diagnostic, so the two are comparable bin-for-bin), the AO -> atom map, the AO-sized
!>        scratch and the trace. Must run after create_rtbse_env (which reads the keyword and
!>        resolves the RI-RS flag) and before initialize_maximum_timestep, whose
!>        reallocate_ft_traces resizes the trace if ENFORCE_MAX_DT changes the step count, and
!>        before initialize_moments, which retains the AO operators only when the flag is set.
!> \param rtbse_env Entry point of the calculation - contains current state of variables
! **************************************************************************************************
   SUBROUTINE build_shell_dipole_maps(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_shell_dipole_maps'

      INTEGER                                            :: handle, iao, iatom, natom, total_sgf
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
            IF (n == 1) THEN
               ! full mantissa: the sum rule against MOMENTS is gated at round-off, and a
               ! 7-significant-digit column would cap the achievable agreement at ~1e-7
               WRITE (print_unit, '(A)') "#      Time [fs]  bin   R [Ang]"// &
                  "                    x,re                    x,im                    y,re"// &
                  "                    y,im                    z,re                    z,im"
               WRITE (print_unit, '(A)') "# the block below is the pre-kick static reference profile"
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

END MODULE rt_bse_shell_dipole
