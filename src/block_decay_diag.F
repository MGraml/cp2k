!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block-decay diagnostic: bins DBCSR block norms by inter-block distance, so the decay
!>        envelope stays visible where the occupation is saturated at 1.0. Own module because
!>        src/emd/ sits above src/gw_* in the dependency stack.
!> \author Max Graml
! **************************************************************************************************
MODULE block_decay_diag

   USE cp_dbcsr_api,                    ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: decay_nrung,&
                                              post_scf_bandstructure_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_decay_diag'

   ! Readers select columns by these names, never by position.
   INTEGER, PARAMETER, PUBLIC :: decay_nmat = 6
   CHARACTER(LEN=16), PARAMETER, PUBLIC :: decay_slot_names(decay_nmat) = &
                                           ['rho_AO          ', 'rho_grid        ', 'rhoW_grid       ', &
                                            'sigma_AO        ', 'W0_grid         ', 'W0_AO           ']

   PUBLIC :: accumulate_block_decay, &
             build_block_decay_maps, &
             print_block_decay, &
             reduce_block_decay, &
             zero_block_decay

CONTAINS

! **************************************************************************************************
!> \brief Size the distance bins and the geometric all-pairs denominator (once). Block positions come
!>        straight from the persisted caches: an AO block is one atom, and a grid block is
!>        atom-aligned, so its centroid is a well-defined position.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE build_block_decay_maps(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_block_decay_maps'

      INTEGER                                            :: handle, i, ibin, j, natom, nbins, nchunks
      REAL(KIND=dp)                                      :: r_ij
      REAL(KIND=dp), DIMENSION(3)                        :: span

      IF (bs_env%ri_rs%block_decay_maps_built) RETURN
      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(bs_env%ri_rs%atom_centers))
      CPASSERT(ALLOCATED(bs_env%ri_rs%chunk_centroids))
      natom = SIZE(bs_env%ri_rs%atom_centers, 2)
      nchunks = SIZE(bs_env%ri_rs%chunk_centroids, 2)

      ! 1 Angstrom bins, in Bohr; enough of them to span the atom bounding-box diagonal.
      IF (bs_env%ri_rs%decay_bin_width <= 0.0_dp) bs_env%ri_rs%decay_bin_width = 1.0_dp/angstrom
      span(:) = MAXVAL(bs_env%ri_rs%atom_centers, DIM=2) - MINVAL(bs_env%ri_rs%atom_centers, DIM=2)
      nbins = INT(SQRT(SUM(span(:)**2))/bs_env%ri_rs%decay_bin_width) + 2
      bs_env%ri_rs%decay_nbins = nbins

      ! source= so the real-input entry points, which never pass through the complex wrapper's
      ! zero_block_decay, still start from a defined state
      ALLOCATE (bs_env%ri_rs%decay_maxnorm(nbins, decay_nmat, 2), source=0.0_dp)
      ALLOCATE (bs_env%ri_rs%decay_count(nbins, decay_nmat, 2), source=0)
      ALLOCATE (bs_env%ri_rs%decay_histcount(nbins, decay_nmat, 2, decay_nrung), source=0)
      ALLOCATE (bs_env%ri_rs%decay_npairs(nbins, 2), source=0)

      ! Pure geometry (no DBCSR): what turns histcount into a fraction of the pairs that exist at R,
      ! rather than "survivors of whatever eps this run happened to use".
      DO i = 1, natom
         DO j = 1, natom
            r_ij = NORM2(bs_env%ri_rs%atom_centers(:, i) - bs_env%ri_rs%atom_centers(:, j))
            ibin = MIN(INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1, nbins)
            bs_env%ri_rs%decay_npairs(ibin, 1) = bs_env%ri_rs%decay_npairs(ibin, 1) + 1
         END DO
      END DO
      DO i = 1, nchunks
         DO j = 1, nchunks
            r_ij = NORM2(bs_env%ri_rs%chunk_centroids(:, i) - bs_env%ri_rs%chunk_centroids(:, j))
            ibin = MIN(INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1, nbins)
            bs_env%ri_rs%decay_npairs(ibin, 2) = bs_env%ri_rs%decay_npairs(ibin, 2) + 1
         END DO
      END DO

      bs_env%ri_rs%block_decay_maps_built = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE build_block_decay_maps

! **************************************************************************************************
!> \brief Bin a matrix's blocks: max block Frobenius norm plus eps-ladder survival counts vs
!>        inter-block distance. Rank-local; the caller reduces.
!> \param matrix ...
!> \param imat slot index, 1 rho_AO, 2 rho_grid, 3 rhoW_grid, 4 sigma_AO, 5 W0_grid, 6 W0_AO
!> \param is_grid .TRUE. for grid x grid (position = chunk centroid), .FALSE. for AO x AO (atom)
!> \param bs_env ...
!> \param ipart 1 = Re (default), 2 = Im
! **************************************************************************************************
   SUBROUTINE accumulate_block_decay(matrix, imat, is_grid, bs_env, ipart)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: imat
      LOGICAL, INTENT(IN)                                :: is_grid
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN), OPTIONAL                      :: ipart

      INTEGER                                            :: col, ibin, irung, ivar, row
      REAL(KIND=dp)                                      :: nrm, r_ij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      ivar = 1
      IF (PRESENT(ipart)) ivar = ipart

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk)
         IF (is_grid) THEN
            r_ij = NORM2(bs_env%ri_rs%chunk_centroids(:, row) - bs_env%ri_rs%chunk_centroids(:, col))
         ELSE
            r_ij = NORM2(bs_env%ri_rs%atom_centers(:, row) - bs_env%ri_rs%atom_centers(:, col))
         END IF
         ibin = MIN(INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1, bs_env%ri_rs%decay_nbins)
         nrm = NORM2(blk)
         bs_env%ri_rs%decay_maxnorm(ibin, imat, ivar) = MAX(bs_env%ri_rs%decay_maxnorm(ibin, imat, ivar), nrm)
         bs_env%ri_rs%decay_count(ibin, imat, ivar) = bs_env%ri_rs%decay_count(ibin, imat, ivar) + 1
         DO irung = 1, decay_nrung
            IF (nrm >= bs_env%ri_rs%decay_eps_ladder(irung)) THEN
               bs_env%ri_rs%decay_histcount(ibin, imat, ivar, irung) = &
                  bs_env%ri_rs%decay_histcount(ibin, imat, ivar, irung) + 1
            END IF
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE accumulate_block_decay

! **************************************************************************************************
!> \brief Zero the histograms. Separate from accumulation because the complex kernel runs Re and Im
!>        passes that must share one set of bins.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE zero_block_decay(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      bs_env%ri_rs%decay_maxnorm(:, :, :) = 0.0_dp
      bs_env%ri_rs%decay_count(:, :, :) = 0
      bs_env%ri_rs%decay_histcount(:, :, :, :) = 0
      ! decay_npairs is geometric, identical on every rank, set once in build_block_decay_maps

   END SUBROUTINE zero_block_decay

! **************************************************************************************************
!> \brief Reduce the histograms across ranks. Collective - every rank must call it, so it cannot
!>        live inside the rank-0 print guard.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE reduce_block_decay(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: ivar

      ! per variant: mp_max has no rank-3 overload; the histogram sums whole-array (rank-4)
      DO ivar = 1, 2
         CALL bs_env%para_env%max(bs_env%ri_rs%decay_maxnorm(:, :, ivar))
         CALL bs_env%para_env%sum(bs_env%ri_rs%decay_count(:, :, ivar))
      END DO
      CALL bs_env%para_env%sum(bs_env%ri_rs%decay_histcount)

   END SUBROUTINE reduce_block_decay

! **************************************************************************************************
!> \brief Emit the per-step decay table: a DECAY| header, the one-shot geometric NPAIRS| denominators,
!>        then DECAY_RE| / DECAY_IM| rows. Rank-0 only; the caller reduces first.
!> \param bs_env ...
!> \param unit_nr ...
!> \param istep ...
! **************************************************************************************************
   SUBROUTINE print_block_decay(bs_env, unit_nr, istep)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: unit_nr, istep

      CHARACTER(LEN=5), PARAMETER :: &
         rung_names(decay_nrung) = ['1e-12', '1e-10', '1e-8 ', '1e-6 ', '1e-4 ', '1e-2 ']

      CHARACTER(LEN=10)                                  :: row_tag
      CHARACTER(LEN=2048)                                :: header
      INTEGER                                            :: ibin, imat, irung, ivar
      REAL(KIND=dp)                                      :: r_ang

      IF (.NOT. bs_env%ri_rs%decay_hdr_printed) THEN
         bs_env%ri_rs%decay_hdr_printed = .TRUE.
         header = ' DECAY|   step  bin  R[Ang]'
         DO imat = 1, decay_nmat
            header = TRIM(header)//' max_'//TRIM(decay_slot_names(imat))
            DO irung = 1, decay_nrung
               header = TRIM(header)//' n>='//TRIM(rung_names(irung))//'_'//TRIM(decay_slot_names(imat))
            END DO
         END DO
         WRITE (unit_nr, '(A)') TRIM(header)//' n_blk'
         CALL print_block_decay_npairs(bs_env, unit_nr)
      END IF

      DO ivar = 1, 2
         ! skip a part that was never accumulated
         IF (SUM(bs_env%ri_rs%decay_count(:, :, ivar)) == 0) CYCLE
         IF (ivar == 1) THEN
            row_tag = 'DECAY_RE| '
         ELSE
            row_tag = 'DECAY_IM| '
         END IF
         DO ibin = 1, bs_env%ri_rs%decay_nbins
            IF (ALL(bs_env%ri_rs%decay_count(ibin, :, ivar) == 0)) CYCLE
            r_ang = (REAL(ibin, dp) - 0.5_dp)*bs_env%ri_rs%decay_bin_width*angstrom
            WRITE (unit_nr, '(1X,A10,I7,I5,F8.2,6(ES14.6,6I12),I12)') &
               row_tag, istep, ibin, r_ang, &
               (bs_env%ri_rs%decay_maxnorm(ibin, imat, ivar), &
                bs_env%ri_rs%decay_histcount(ibin, imat, ivar, 1:decay_nrung), imat=1, decay_nmat), &
               MAXVAL(bs_env%ri_rs%decay_count(ibin, :, ivar))
         END DO
      END DO

   END SUBROUTINE print_block_decay

! **************************************************************************************************
!> \brief Emit the geometric n_pairs tables once: the occupancy denominators, rank-identical data.
!> \param bs_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_block_decay_npairs(bs_env, unit_nr)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: ibin
      REAL(KIND=dp)                                      :: r_ang

      IF (bs_env%ri_rs%decay_npairs_printed) RETURN
      bs_env%ri_rs%decay_npairs_printed = .TRUE.

      WRITE (unit_nr, '(A)') ' NPAIRS|  bin  R[Ang]     n_pairs'
      DO ibin = 1, bs_env%ri_rs%decay_nbins
         IF (bs_env%ri_rs%decay_npairs(ibin, 1) /= 0) THEN
            r_ang = (REAL(ibin, dp) - 0.5_dp)*bs_env%ri_rs%decay_bin_width*angstrom
            WRITE (unit_nr, '(1X,A13,I5,F8.2,I12)') 'NPAIRS_ATOM| ', ibin, r_ang, bs_env%ri_rs%decay_npairs(ibin, 1)
         END IF
      END DO
      DO ibin = 1, bs_env%ri_rs%decay_nbins
         IF (bs_env%ri_rs%decay_npairs(ibin, 2) /= 0) THEN
            r_ang = (REAL(ibin, dp) - 0.5_dp)*bs_env%ri_rs%decay_bin_width*angstrom
            WRITE (unit_nr, '(1X,A13,I5,F8.2,I12)') 'NPAIRS_CHUNK|', ibin, r_ang, bs_env%ri_rs%decay_npairs(ibin, 2)
         END IF
      END DO

   END SUBROUTINE print_block_decay_npairs

END MODULE block_decay_diag
