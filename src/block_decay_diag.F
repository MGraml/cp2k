!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block-decay diagnostic: bins DBCSR block norms by inter-atomic distance so the decay
!>        envelope stays visible where the occupation is saturated at 1.0. Shared by the RI-RS
!>        GW self-energy and the linearized RT-BSE kernel, hence its own module: src/emd/ sits
!>        above src/gw_* in the dependency stack.
!> \author Max Graml
! **************************************************************************************************
MODULE block_decay_diag

   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_iterator_blocks_left,&
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

   ! Descriptive per-slot column names for the emitted decay tables; readers parse by header name.
   ! Slots 1-8 RT-BSE per-step/static, 9-20 GW-side capture, 21-23 RT-BSE static controls.
   CHARACTER(LEN=16), PARAMETER, PUBLIC :: decay_slot_names(23) = &
                                           ['rho_AO          ', 'rho_grid        ', 'rhoW_grid       ', 'sigma_AO        ', &
                                            'W0_grid         ', 'W0_AO           ', 'W_PQ            ', 'D               ', &
                                            'G_occ_AO        ', 'G_vir_AO        ', 'G_occ_grid      ', 'G_vir_grid      ', &
                                            'W_grid          ', 'Sigma_c_neg_grid', 'Sigma_c_pos_grid', 'Sigma_c_neg_AO  ', &
                                            'Sigma_x_grid    ', 'Sigma_x_AO      ', 'Sigma_c_pos_AO  ', 'W_AO            ', &
                                            'Dgrid           ', 'Gvir_AO         ', 'Gvir_grid       ']

   PUBLIC :: accumulate_block_decay, &
             build_block_decay_maps, &
             decay_header_append_slot, &
             print_block_decay

CONTAINS

! **************************************************************************************************
!> \brief Build the chunk->atom maps for the block-decay diagnostic (once). An AO block indexes its
!>        atom directly; a grid chunk (<=32 kind-major-ordered grid points) maps to the atom owning its
!>        midpoint, flagged straddling when its first and last points fall in different atoms. Also
!>        sizes the 1 A distance bins from the atom bounding-box diagonal.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE build_block_decay_maps(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_block_decay_maps'

      INTEGER                                            :: cend, cmid, cstart, handle, iatom, ibin, &
                                                            ichunk, ipt, jatom, jchunk, natom, &
                                                            nbins, nchunks, total_npts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: chunk_start, pt_to_atom
      INTEGER, DIMENSION(:), POINTER                     :: blk_grid_sizes
      REAL(KIND=dp)                                      :: diag, r_ij
      REAL(KIND=dp), DIMENSION(3)                        :: span

      IF (bs_env%ri_rs%block_decay_maps_built) RETURN
      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, row_blk_size=blk_grid_sizes)
      nchunks = SIZE(blk_grid_sizes)

      ! first grid-point index of each chunk
      ALLOCATE (chunk_start(nchunks + 1))
      chunk_start(1) = 1
      DO ichunk = 1, nchunks
         chunk_start(ichunk + 1) = chunk_start(ichunk) + blk_grid_sizes(ichunk)
      END DO

      ! Grid points are laid out kind-major, so grid_offsets does NOT ascend with the atom index
      ! whenever kinds interleave in the coordinate file. Mark each atom's start and forward-fill,
      ! which assumes no ordering at all.
      natom = SIZE(bs_env%ri_rs%grid_offsets) - 1
      total_npts = bs_env%ri_rs%grid_offsets(natom + 1) - 1
      CPASSERT(chunk_start(nchunks + 1) - 1 == total_npts)

      ALLOCATE (pt_to_atom(total_npts))
      pt_to_atom(:) = 0
      DO iatom = 1, natom
         pt_to_atom(bs_env%ri_rs%grid_offsets(iatom)) = iatom
      END DO
      iatom = 0
      DO ipt = 1, total_npts
         IF (pt_to_atom(ipt) > 0) iatom = pt_to_atom(ipt)
         pt_to_atom(ipt) = iatom
      END DO
      CPASSERT(ALL(pt_to_atom > 0))

      ALLOCATE (bs_env%ri_rs%chunk_to_atom(nchunks), bs_env%ri_rs%chunk_straddles(nchunks))
      DO ichunk = 1, nchunks
         cstart = chunk_start(ichunk)
         cend = chunk_start(ichunk + 1) - 1
         cmid = (cstart + cend)/2
         bs_env%ri_rs%chunk_to_atom(ichunk) = pt_to_atom(cmid)
         bs_env%ri_rs%chunk_straddles(ichunk) = (pt_to_atom(cstart) /= pt_to_atom(cend))
      END DO
      DEALLOCATE (pt_to_atom, chunk_start)

      span(:) = MAXVAL(bs_env%ri_rs%atom_pos, DIM=2) - MINVAL(bs_env%ri_rs%atom_pos, DIM=2)
      diag = SQRT(SUM(span(:)**2))
      nbins = INT(diag/bs_env%ri_rs%decay_bin_width) + 2
      bs_env%ri_rs%decay_nbins = nbins
      ALLOCATE (bs_env%ri_rs%decay_maxnorm(nbins, bs_env%ri_rs%decay_nmat, 2), &
                bs_env%ri_rs%decay_count(nbins, bs_env%ri_rs%decay_nmat, 2), &
                bs_env%ri_rs%decay_histcount(nbins, bs_env%ri_rs%decay_nmat, 2, decay_nrung), &
                bs_env%ri_rs%decay_npairs(nbins, 2))

      ! Geometric all-pairs denominator n_pairs(R), pure geometry (no DBCSR): what makes
      ! occupancy(R; eps) = histcount/n_pairs a fraction rather than "survivors of this run's eps".
      bs_env%ri_rs%decay_npairs(:, :) = 0
      DO iatom = 1, natom
         DO jatom = 1, natom
            r_ij = NORM2(bs_env%ri_rs%atom_pos(:, iatom) - bs_env%ri_rs%atom_pos(:, jatom))
            ibin = INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1
            IF (ibin > nbins) ibin = nbins
            bs_env%ri_rs%decay_npairs(ibin, 1) = bs_env%ri_rs%decay_npairs(ibin, 1) + 1
         END DO
      END DO
      DO ichunk = 1, nchunks
         DO jchunk = 1, nchunks
            r_ij = NORM2(bs_env%ri_rs%atom_pos(:, bs_env%ri_rs%chunk_to_atom(ichunk)) - &
                         bs_env%ri_rs%atom_pos(:, bs_env%ri_rs%chunk_to_atom(jchunk)))
            ibin = INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1
            IF (ibin > nbins) ibin = nbins
            bs_env%ri_rs%decay_npairs(ibin, 2) = bs_env%ri_rs%decay_npairs(ibin, 2) + 1
         END DO
      END DO

      bs_env%ri_rs%block_decay_maps_built = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE build_block_decay_maps

! **************************************************************************************************
!> \brief Bin a matrix's blocks into the decay histogram: max block Frobenius norm plus eps-ladder
!>        survival counts vs inter-atomic distance. Every block feeds the ALL variant (straddlers
!>        at their majority atom); only non-straddling blocks feed the CLEAN variant. Rank-local;
!>        the caller reduces.
!> \param matrix ...
!> \param imat matched-pair index (1 rho_ao, 2 rho_grid, 3 rhoW_grid, 4 sigma_ao, 5 W0_grid, 6 W0_ao)
!> \param is_grid .TRUE. for grid×grid (block index = chunk), .FALSE. for AO×AO (block index = atom)
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE accumulate_block_decay(matrix, imat, is_grid, bs_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: imat
      LOGICAL, INTENT(IN)                                :: is_grid
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: atom_i, atom_j, col, ibin, irung, row
      LOGICAL                                            :: straddles
      REAL(KIND=dp)                                      :: nrm, r_ij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk)
         IF (is_grid) THEN
            atom_i = bs_env%ri_rs%chunk_to_atom(row)
            atom_j = bs_env%ri_rs%chunk_to_atom(col)
            straddles = bs_env%ri_rs%chunk_straddles(row) .OR. bs_env%ri_rs%chunk_straddles(col)
         ELSE
            atom_i = row
            atom_j = col
            straddles = .FALSE.
         END IF
         r_ij = NORM2(bs_env%ri_rs%atom_pos(:, atom_i) - bs_env%ri_rs%atom_pos(:, atom_j))
         ibin = INT(r_ij/bs_env%ri_rs%decay_bin_width) + 1
         IF (ibin > bs_env%ri_rs%decay_nbins) ibin = bs_env%ri_rs%decay_nbins
         nrm = NORM2(blk)
         bs_env%ri_rs%decay_maxnorm(ibin, imat, 1) = MAX(bs_env%ri_rs%decay_maxnorm(ibin, imat, 1), nrm)
         bs_env%ri_rs%decay_count(ibin, imat, 1) = bs_env%ri_rs%decay_count(ibin, imat, 1) + 1
         DO irung = 1, decay_nrung
            IF (nrm >= bs_env%ri_rs%decay_eps_ladder(irung)) THEN
               bs_env%ri_rs%decay_histcount(ibin, imat, 1, irung) = &
                  bs_env%ri_rs%decay_histcount(ibin, imat, 1, irung) + 1
            END IF
         END DO
         IF (.NOT. straddles) THEN
            bs_env%ri_rs%decay_maxnorm(ibin, imat, 2) = MAX(bs_env%ri_rs%decay_maxnorm(ibin, imat, 2), nrm)
            bs_env%ri_rs%decay_count(ibin, imat, 2) = bs_env%ri_rs%decay_count(ibin, imat, 2) + 1
            DO irung = 1, decay_nrung
               IF (nrm >= bs_env%ri_rs%decay_eps_ladder(irung)) THEN
                  bs_env%ri_rs%decay_histcount(ibin, imat, 2, irung) = &
                     bs_env%ri_rs%decay_histcount(ibin, imat, 2, irung) + 1
               END IF
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE accumulate_block_decay

! **************************************************************************************************
!> \brief Append one slot's column names (max_<name> then n>=<eps>_<name> per ladder rung) to a
!>        header line under construction; readers select columns by these names, never by position.
!> \param header ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE decay_header_append_slot(header, name)
      CHARACTER(LEN=*), INTENT(INOUT)                    :: header
      CHARACTER(LEN=*), INTENT(IN)                       :: name

      CHARACTER(LEN=5), PARAMETER :: &
         rung_names(decay_nrung) = ['1e-12', '1e-10', '1e-8 ', '1e-6 ', '1e-4 ', '1e-2 ']

      INTEGER                                            :: irung

      header = TRIM(header)//' max_'//TRIM(name)
      DO irung = 1, decay_nrung
         header = TRIM(header)//' n>='//TRIM(rung_names(irung))//'_'//TRIM(name)
      END DO

   END SUBROUTINE decay_header_append_slot

! **************************************************************************************************
!> \brief Emit the per-step block-decay tables. See print_block_decay_static (rt_bse_ri_rs) for the
!>        step-independent static table.
!> \param bs_env ...
!> \param unit_nr ...
!> \param istep ...
! **************************************************************************************************
   SUBROUTINE print_block_decay(bs_env, unit_nr, istep)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: unit_nr, istep

      CHARACTER(LEN=10)                                  :: tag
      CHARACTER(LEN=1024)                                :: header
      INTEGER                                            :: ibin, imat, ivar
      REAL(KIND=dp)                                      :: r_ang

      IF (istep == 0) THEN
         header = ' DECAY|   step  bin  R[Ang]'
         DO imat = 1, 6
            CALL decay_header_append_slot(header, decay_slot_names(imat))
         END DO
         header = TRIM(header)//' n_all'
         WRITE (unit_nr, '(A)') TRIM(header)
      END IF
      DO ivar = 1, 2
         IF (ivar == 1) THEN
            tag = 'DECAY_ALL|'
         ELSE
            tag = 'DECAY_CLN|'
         END IF
         DO ibin = 1, bs_env%ri_rs%decay_nbins
            IF (bs_env%ri_rs%decay_count(ibin, 1, ivar) == 0 .AND. &
                bs_env%ri_rs%decay_count(ibin, 2, ivar) == 0) CYCLE
            r_ang = (REAL(ibin, dp) - 0.5_dp)*bs_env%ri_rs%decay_bin_width*angstrom
            WRITE (unit_nr, '(1X,A10,I7,I5,F8.2,6(ES14.6,6I12),I12)') &
               tag, istep, ibin, r_ang, &
               (bs_env%ri_rs%decay_maxnorm(ibin, imat, ivar), &
                bs_env%ri_rs%decay_histcount(ibin, imat, ivar, 1:decay_nrung), imat=1, 6), &
               bs_env%ri_rs%decay_count(ibin, 1, ivar)
         END DO
      END DO

   END SUBROUTINE print_block_decay

END MODULE block_decay_diag
