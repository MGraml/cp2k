!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   DBCSR output in CP2K
!> \author  VW
!> \date    2009-09-09
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2009-09-09
! **************************************************************************************************
MODULE cp_dbcsr_output
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_get_data_size, dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: symmetrize_matrix
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: nso
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_output'

   PUBLIC :: cp_dbcsr_write_sparse_matrix
   PUBLIC :: cp_dbcsr_write_matrix_dist
   PUBLIC :: write_fm_with_basis_info

   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief Print a spherical matrix of blacs type.
!> \param blacs_matrix ...
!> \param before ...
!> \param after ...
!> \param qs_env ...
!> \param para_env ...
!> \param first_row ...
!> \param last_row ...
!> \param first_col ...
!> \param last_col ...
!> \param output_unit ...
!> \param omit_headers Write only the matrix data, not the row/column headers
!> \author Creation (12.06.2001,MK)
!>       Allow for printing of a sub-matrix (01.07.2003,MK)
! **************************************************************************************************
   SUBROUTINE write_fm_with_basis_info(blacs_matrix, before, after, qs_env, para_env, &
                                       first_row, last_row, first_col, last_col, output_unit, omit_headers)

      TYPE(cp_fm_type), INTENT(IN)                       :: blacs_matrix
      INTEGER, INTENT(IN)                                :: before, after
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: first_row, last_row, first_col, last_col
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN), OPTIONAL                      :: omit_headers

      CHARACTER(LEN=60)                                  :: matrix_name
      INTEGER                                            :: col1, col2, ncol_global, nrow_global, &
                                                            nsgf, row1, row2
      LOGICAL                                            :: my_omit_headers
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      IF (.NOT. ASSOCIATED(blacs_matrix%matrix_struct)) RETURN
      CALL cp_fm_get_info(blacs_matrix, name=matrix_name, nrow_global=nrow_global, &
                          ncol_global=ncol_global)

      ALLOCATE (matrix(nrow_global, ncol_global))
      CALL cp_fm_get_submatrix(blacs_matrix, matrix)

      ! *** Get the matrix dimension and check the optional arguments ***
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, nsgf=nsgf)

      IF (PRESENT(first_row)) THEN
         row1 = MAX(1, first_row)
      ELSE
         row1 = 1
      END IF

      IF (PRESENT(last_row)) THEN
         row2 = MIN(nsgf, last_row)
      ELSE
         row2 = nsgf
      END IF

      IF (PRESENT(first_col)) THEN
         col1 = MAX(1, first_col)
      ELSE
         col1 = 1
      END IF

      IF (PRESENT(last_col)) THEN
         col2 = MIN(nsgf, last_col)
      ELSE
         col2 = nsgf
      END IF

      IF (PRESENT(omit_headers)) THEN
         my_omit_headers = omit_headers
      ELSE
         my_omit_headers = .FALSE.
      END IF

      CALL write_matrix_sym(matrix, matrix_name, before, after, qs_env, para_env, &
                            row1, row2, col1, col2, output_unit, omit_headers=my_omit_headers)

      ! *** Release work storage ***
      IF (ASSOCIATED(matrix)) THEN
         DEALLOCATE (matrix)
      END IF

   END SUBROUTINE write_fm_with_basis_info

! **************************************************************************************************
!> \brief ...
!> \param sparse_matrix ...
!> \param before ...
!> \param after ...
!> \param qs_env ...
!> \param para_env ...
!> \param first_row ...
!> \param last_row ...
!> \param first_col ...
!> \param last_col ...
!> \param scale ...
!> \param output_unit ...
!> \param omit_headers Write only the matrix data, not the row/column headers
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_write_sparse_matrix(sparse_matrix, before, after, qs_env, para_env, &
                                           first_row, last_row, first_col, last_col, scale, &
                                           output_unit, omit_headers)

      TYPE(dbcsr_type)                                   :: sparse_matrix
      INTEGER, INTENT(IN)                                :: before, after
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: first_row, last_row, first_col, last_col
      REAL(dp), INTENT(IN), OPTIONAL                     :: scale
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN), OPTIONAL                      :: omit_headers

      CHARACTER(LEN=default_string_length)               :: matrix_name
      INTEGER                                            :: col1, col2, dim_col, dim_row, row1, row2
      LOGICAL                                            :: my_omit_headers, print_sym
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (matrix)

      CALL copy_repl_dbcsr_to_repl_fm(sparse_matrix, matrix)

      CALL para_env%sum(matrix)

      SELECT CASE (dbcsr_get_matrix_type(sparse_matrix))
      CASE (dbcsr_type_symmetric)
         CALL symmetrize_matrix(matrix, "upper_to_lower")
         print_sym = .TRUE.
      CASE (dbcsr_type_antisymmetric)
         CALL symmetrize_matrix(matrix, "anti_upper_to_lower")
         print_sym = .TRUE.
      CASE (dbcsr_type_no_symmetry)
         print_sym = .FALSE.
      CASE DEFAULT
         CPABORT("WRONG")
      END SELECT

      ! *** Get the matrix dimension and check the optional arguments ***
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      dim_row = SIZE(matrix, 1)
      dim_col = SIZE(matrix, 2)

      IF (PRESENT(first_row)) THEN
         row1 = MAX(1, first_row)
      ELSE
         row1 = 1
      END IF

      IF (PRESENT(last_row)) THEN
         row2 = MIN(dim_row, last_row)
      ELSE
         row2 = dim_row
      END IF

      IF (PRESENT(first_col)) THEN
         col1 = MAX(1, first_col)
      ELSE
         col1 = 1
      END IF

      IF (PRESENT(last_col)) THEN
         col2 = MIN(dim_col, last_col)
      ELSE
         col2 = dim_col
      END IF

      IF (PRESENT(scale)) THEN
         matrix = matrix*scale
      END IF

      IF (PRESENT(omit_headers)) THEN
         my_omit_headers = omit_headers
      ELSE
         my_omit_headers = .FALSE.
      END IF

      CALL dbcsr_get_info(sparse_matrix, name=matrix_name)
      IF (print_sym) THEN
         CALL write_matrix_sym(matrix, matrix_name, before, after, qs_env, para_env, &
                               row1, row2, col1, col2, output_unit, my_omit_headers)
      ELSE
         CALL write_matrix_gen(matrix, matrix_name, before, after, para_env, &
                               row1, row2, col1, col2, output_unit, my_omit_headers)
      END IF

      IF (ASSOCIATED(matrix)) THEN
         DEALLOCATE (matrix)
      END IF

   END SUBROUTINE cp_dbcsr_write_sparse_matrix

! **************************************************************************************************
!> \brief ...
!> \param sparse_matrix ...
!> \param fm ...
! **************************************************************************************************
   SUBROUTINE copy_repl_dbcsr_to_repl_fm(sparse_matrix, fm)

      TYPE(dbcsr_type)                                   :: sparse_matrix
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fm

      CHARACTER(len=*), PARAMETER :: routineN = 'copy_repl_dbcsr_to_repl_fm'

      INTEGER                                            :: col, handle, i, j, nblkcols_total, &
                                                            nblkrows_total, nc, nr, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: c_offset, r_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(fm)) DEALLOCATE (fm)

      CALL dbcsr_get_info(matrix=sparse_matrix, &
                          col_blk_size=col_blk_size, &
                          row_blk_size=row_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total)

      !> this should be precomputed somewhere else
      ALLOCATE (r_offset(nblkrows_total), c_offset(nblkcols_total))

      r_offset(1) = 1
      DO row = 2, nblkrows_total
         r_offset(row) = r_offset(row - 1) + row_blk_size(row - 1)
      END DO
      nr = SUM(row_blk_size)
      c_offset(1) = 1
      DO col = 2, nblkcols_total
         c_offset(col) = c_offset(col - 1) + col_blk_size(col - 1)
      END DO
      nc = SUM(col_blk_size)
      !<

      ALLOCATE (fm(nr, nc))

      fm(:, :) = 0.0_dp

      CALL dbcsr_iterator_start(iter, sparse_matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         DO j = 1, SIZE(block, 2)
         DO i = 1, SIZE(block, 1)
            fm(r_offset(row) + i - 1, c_offset(col) + j - 1) = block(i, j)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (r_offset, c_offset)

      CALL timestop(handle)

   END SUBROUTINE copy_repl_dbcsr_to_repl_fm

! **************************************************************************************************
!> \brief Write a matrix or a sub-matrix to the output unit (symmetric)
!> \param matrix ...
!> \param matrix_name ...
!> \param before ...
!> \param after ...
!> \param qs_env ...
!> \param para_env ...
!> \param first_row ...
!> \param last_row ...
!> \param first_col ...
!> \param last_col ...
!> \param output_unit ...
!> \param omit_headers Write only the matrix data, not the row/column headers
!> \author Creation (01.07.2003,MK)
! **************************************************************************************************
   SUBROUTINE write_matrix_sym(matrix, matrix_name, before, after, qs_env, para_env, &
                               first_row, last_row, first_col, last_col, output_unit, omit_headers)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name
      INTEGER, INTENT(IN)                                :: before, after
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: first_row, last_row, first_col, &
                                                            last_col, output_unit
      LOGICAL, INTENT(IN)                                :: omit_headers

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=25)                                  :: fmtstr1
      CHARACTER(LEN=35)                                  :: fmtstr2
      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: sgf_symbol
      INTEGER                                            :: from, iatom, icol, ikind, irow, iset, &
                                                            isgf, ishell, iso, jcol, l, left, &
                                                            natom, ncol, ndigits, nset, nsgf, &
                                                            right, to, width
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      IF (output_unit > 0) THEN
         CALL m_flush(output_unit)

         CALL get_qs_env(qs_env=qs_env, &
                         qs_kind_set=qs_kind_set, &
                         atomic_kind_set=atomic_kind_set, &
                         particle_set=particle_set)

         natom = SIZE(particle_set)

         CALL get_qs_kind_set(qs_kind_set=qs_kind_set, nsgf=nsgf)

         ALLOCATE (first_sgf(natom))
         ALLOCATE (last_sgf(natom))
         CALL get_particle_set(particle_set, qs_kind_set, &
                               first_sgf=first_sgf, &
                               last_sgf=last_sgf)

         ! *** Definition of the variable formats ***
         fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
         IF (omit_headers) THEN
            fmtstr2 = "(T2,   (1X,F  .  ))"
         ELSE
            fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"
         END IF

         ! *** Write headline ***
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(matrix_name)

         ! *** Write the variable format strings ***
         ndigits = after

         width = before + ndigits + 3
         ncol = INT(56/width)

         right = MAX((ndigits - 2), 1)
         left = width - right - 5

         WRITE (UNIT=fmtstr1(11:12), FMT="(I2)") ncol
         WRITE (UNIT=fmtstr1(14:15), FMT="(I2)") left
         WRITE (UNIT=fmtstr1(21:22), FMT="(I2)") right

         IF (omit_headers) THEN
            WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") ncol
            WRITE (UNIT=fmtstr2(13:14), FMT="(I2)") width - 1
            WRITE (UNIT=fmtstr2(16:17), FMT="(I2)") ndigits
         ELSE
            WRITE (UNIT=fmtstr2(22:23), FMT="(I2)") ncol
            WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") width - 1
            WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
         END IF

         ! *** Write the matrix in the selected format ***
         DO icol = first_col, last_col, ncol
            from = icol
            to = MIN((from + ncol - 1), last_col)
            IF (.NOT. omit_headers) THEN
               WRITE (UNIT=output_unit, FMT=fmtstr1) (jcol, jcol=from, to)
            END IF
            irow = 1
            DO iatom = 1, natom
               NULLIFY (orb_basis_set)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, &
                                    kind_number=ikind, element_symbol=element_symbol)
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
               IF (ASSOCIATED(orb_basis_set)) THEN
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         nset=nset, nshell=nshell, l=lshell, sgf_symbol=sgf_symbol)
                  isgf = 1
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        l = lshell(ishell, iset)
                        DO iso = 1, nso(l)
                           IF ((irow >= first_row) .AND. (irow <= last_row)) THEN
                              IF (omit_headers) THEN
                                 WRITE (UNIT=output_unit, FMT=fmtstr2) &
                                    (matrix(irow, jcol), jcol=from, to)
                              ELSE
                                 WRITE (UNIT=output_unit, FMT=fmtstr2) &
                                    irow, iatom, element_symbol, sgf_symbol(isgf), &
                                    (matrix(irow, jcol), jcol=from, to)
                              END IF
                           END IF
                           isgf = isgf + 1
                           irow = irow + 1
                        END DO
                     END DO
                  END DO
                  IF ((irow >= first_row) .AND. (irow <= last_row)) THEN
                     WRITE (UNIT=output_unit, FMT="(A)")
                  END IF
               ELSE
                  DO iso = first_sgf(iatom), last_sgf(iatom)
                     IF ((irow >= first_row) .AND. (irow <= last_row)) THEN
                        IF (omit_headers) THEN
                           WRITE (UNIT=output_unit, FMT=fmtstr2) &
                              (matrix(irow, jcol), jcol=from, to)
                        ELSE
                           WRITE (UNIT=output_unit, FMT=fmtstr2) &
                              irow, iatom, element_symbol, " ", &
                              (matrix(irow, jcol), jcol=from, to)
                        END IF
                     END IF
                     irow = irow + 1
                  END DO
                  IF ((irow >= first_row) .AND. (irow <= last_row)) THEN
                     WRITE (UNIT=output_unit, FMT="(A)")
                  END IF
               END IF
            END DO
         END DO

         WRITE (UNIT=output_unit, FMT="(/)")
         DEALLOCATE (first_sgf)
         DEALLOCATE (last_sgf)
      END IF

      CALL para_env%sync()
      IF (output_unit > 0) CALL m_flush(output_unit)

   END SUBROUTINE write_matrix_sym

! **************************************************************************************************
!> \brief Write a matrix not necessarily symmetric (no index with atomic labels)
!> \param matrix ...
!> \param matrix_name ...
!> \param before ...
!> \param after ...
!> \param para_env ...
!> \param first_row ...
!> \param last_row ...
!> \param first_col ...
!> \param last_col ...
!> \param output_unit ...
!> \param omit_headers Write only the matrix data, not the row/column headers
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! **************************************************************************************************
   SUBROUTINE write_matrix_gen(matrix, matrix_name, before, after, para_env, &
                               first_row, last_row, first_col, last_col, output_unit, omit_headers)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name
      INTEGER, INTENT(IN)                                :: before, after
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: first_row, last_row, first_col, &
                                                            last_col, output_unit
      LOGICAL, INTENT(IN)                                :: omit_headers

      CHARACTER(LEN=25)                                  :: fmtstr1
      CHARACTER(LEN=35)                                  :: fmtstr2
      INTEGER                                            :: from, icol, irow, jcol, left, ncol, &
                                                            ndigits, right, to, width

      IF (output_unit > 0) THEN
         CALL m_flush(output_unit)

         ! *** Definition of the variable formats ***
         fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
         IF (omit_headers) THEN
            fmtstr2 = "(T2,   (1X,F  .  ))"
         ELSE
            fmtstr2 = "(T2, I5,        18X,   (1X,F  .  ))"
         END IF

         ! *** Write headline ***
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(matrix_name)

         ! *** Write the variable format strings ***
         ndigits = after

         width = before + ndigits + 3
         ncol = INT(56/width)

         right = MAX((ndigits - 2), 1)
         left = width - right - 5

         WRITE (UNIT=fmtstr1(11:12), FMT="(I2)") ncol
         WRITE (UNIT=fmtstr1(14:15), FMT="(I2)") left
         WRITE (UNIT=fmtstr1(21:22), FMT="(I2)") right

         IF (omit_headers) THEN
            WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") ncol
            WRITE (UNIT=fmtstr2(13:14), FMT="(I2)") width - 1
            WRITE (UNIT=fmtstr2(16:17), FMT="(I2)") ndigits
         ELSE
            WRITE (UNIT=fmtstr2(22:23), FMT="(I2)") ncol
            WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") width - 1
            WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
         END IF

         ! *** Write the matrix in the selected format ***
         DO icol = first_col, last_col, ncol
            from = icol
            to = MIN((from + ncol - 1), last_col)
            IF (.NOT. omit_headers) THEN
               WRITE (UNIT=output_unit, FMT=fmtstr1) (jcol, jcol=from, to)
            END IF
            irow = 1
            DO irow = first_row, last_row
               IF (omit_headers) THEN
                  WRITE (UNIT=output_unit, FMT=fmtstr2) &
                     irow, (matrix(irow, jcol), jcol=from, to)
               ELSE
                  WRITE (UNIT=output_unit, FMT=fmtstr2) &
                     (matrix(irow, jcol), jcol=from, to)
               END IF
            END DO
         END DO

         WRITE (UNIT=output_unit, FMT="(/)")
      END IF

      CALL para_env%sync()
      IF (output_unit > 0) CALL m_flush(output_unit)

   END SUBROUTINE write_matrix_gen

! **************************************************************************************************
!> \brief Print the distribution of a sparse matrix.
!> \param matrix ...
!> \param output_unit ...
!> \param para_env ...
!> \par History
!>      Creation (25.06.2003,MK)
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_write_matrix_dist(matrix, output_unit, para_env)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_write_matrix_dist'
      LOGICAL, PARAMETER                                 :: full_output = .FALSE.

      CHARACTER                                          :: matrix_type
      CHARACTER(LEN=default_string_length)               :: matrix_name
      INTEGER                                            :: handle, ipe, mype, natom, nblock_max, &
                                                            nelement_max, npe, nrow, tmp(2)
      INTEGER(KIND=int_8)                                :: nblock_sum, nblock_tot, nelement_sum
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nblock, nelement
      LOGICAL                                            :: ionode
      REAL(KIND=dp)                                      :: occupation
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      CALL timeset(routineN, handle)

      ionode = para_env%is_source()
      mype = para_env%mepos + 1
      npe = para_env%num_pe

      ! *** Allocate work storage ***
      ALLOCATE (nblock(npe))
      nblock(:) = 0

      ALLOCATE (nelement(npe))
      nelement(:) = 0

      nblock(mype) = dbcsr_get_num_blocks(matrix)
      nelement(mype) = dbcsr_get_data_size(matrix)

      CALL dbcsr_get_info(matrix=matrix, &
                          name=matrix_name, &
                          matrix_type=matrix_type, &
                          nblkrows_total=natom, &
                          nfullrows_total=nrow)

      IF (full_output) THEN
         ! XXXXXXXX should gather/scatter this on ionode
         CALL para_env%sum(nblock)
         CALL para_env%sum(nelement)

         nblock_sum = SUM(INT(nblock, KIND=int_8))
         nelement_sum = SUM(INT(nelement, KIND=int_8))
      ELSE
         nblock_sum = nblock(mype)
         nblock_max = nblock(mype)
         nelement_sum = nelement(mype)
         nelement_max = nelement(mype)
         CALL para_env%sum(nblock_sum)
         CALL para_env%sum(nelement_sum)
         tmp = (/nblock_max, nelement_max/)
         CALL para_env%max(tmp)
         nblock_max = tmp(1); nelement_max = tmp(2)
      END IF

      IF (matrix_type == dbcsr_type_symmetric .OR. &
          matrix_type == dbcsr_type_antisymmetric) THEN
         nblock_tot = INT(natom, KIND=int_8)*INT(natom + 1, KIND=int_8)/2
      ELSE
         nblock_tot = INT(natom, KIND=int_8)**2
      END IF

      occupation = -1.0_dp
      IF (nblock_tot .NE. 0) occupation = 100.0_dp*REAL(nblock_sum, dp)/REAL(nblock_tot, dp)

      IF (ionode) THEN
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "DISTRIBUTION OF THE "//TRIM(matrix_name)
         IF (full_output) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/,/,(I9,T27,I10,T55,I10))") &
               "Process    Number of matrix blocks   Number of matrix elements", &
               (ipe - 1, nblock(ipe), nelement(ipe), ipe=1, npe)
            WRITE (UNIT=output_unit, FMT="(/,T7,A3,T27,I10,T55,I10)") &
               "Sum", nblock_sum, nelement_sum
            WRITE (UNIT=output_unit, FMT="(/,T7,A3,T27,I10,A,F5.1,A,T55,I10,A,F5.1,A)") &
               " of", nblock_tot, " (", occupation, " % occupation)"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T15,A,T68,I13)") "Number  of non-zero blocks:", nblock_sum
            WRITE (UNIT=output_unit, FMT="(T15,A,T75,F6.2)") "Percentage non-zero blocks:", occupation
            WRITE (UNIT=output_unit, FMT="(T15,A,T68,I13)") "Average number of blocks per CPU:", &
               (nblock_sum + npe - 1)/npe
            WRITE (UNIT=output_unit, FMT="(T15,A,T68,I13)") "Maximum number of blocks per CPU:", nblock_max
            WRITE (UNIT=output_unit, FMT="(T15,A,T68,I13)") "Average number of matrix elements per CPU:", &
               (nelement_sum + npe - 1)/npe
            WRITE (UNIT=output_unit, FMT="(T15,A,T68,I13)") "Maximum number of matrix elements per CPU:", &
               nelement_max
         END IF
      END IF

      ! *** Release work storage ***
      DEALLOCATE (nblock)

      DEALLOCATE (nelement)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_write_matrix_dist

END MODULE cp_dbcsr_output
