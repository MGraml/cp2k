!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2017 created [Jan Wilhelm]
! **************************************************************************************************
MODULE bse
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE group_dist_types,                ONLY: get_group_dist,&
                                              group_dist_d1_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse'

   PUBLIC :: mult_B_with_W_and_fill_local_3c_arrays, do_subspace_iterations, deallocate_matrices_BSE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param num_davidson_iter ...
!> \param eps_res ...
!> \param Eigenval ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo, virtual, num_Z_vectors_init, &
                                     num_davidson_iter, eps_res, Eigenval, para_env, &
                                     num_exc_en, num_add_start_z_space, fac_max_z_space, num_new_t)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER                                            :: homo, virtual, num_Z_vectors_init, num_davidson_iter, &
                                                            num_exc_en, num_add_start_z_space, fac_max_z_space, &
                                                            num_new_t, num_Z_vectors
      REAL(KIND=dp)                                      :: eps_res
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'do_subspace_iterations'

      INTEGER                                            :: handle, i_iter, local_RI_size,i_print,j_print,k_print,num_max_z_space
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_ia_tmp, M_ji_tmp, RI_vector, Subspace_eigenvec,AZ_reshaped, Z_vectors_reshaped, &
                                                            Subspace_residuals_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, BZ, Z_vectors, Preconditioner, Subspace_add_dir, W_vectors


      CALL timeset(routineN, handle)
      
      !MG to del
      print *,"num_exc_en",num_exc_en
      print*,"num_add_start_z_space",num_add_start_z_space
      print *,"num_Z_vectors_init", num_Z_vectors_init
      print *,"fac_max_z_space", fac_max_z_space
      print *, "num_new_t", num_new_t
      print *, "eps_res", eps_res
      print *, "num_davidson_iter", num_davidson_iter
      !MG to del: avoid unused errors
      i_iter = 0
      num_Z_vectors = num_Z_vectors_init
      num_max_z_space = num_Z_vectors_init * fac_max_z_space
      ! JW hack 2del
      eps_res = 0.01_dp

      

      !print *, "Printing B_bar_iaQ_bse_local of shape", SHAPE(B_bar_iaQ_bse_local)
      !do i_print=1,SIZE(B_bar_iaQ_bse_local,2)
      !   print *, B_bar_iaQ_bse_local(:,:,1)
      !end do
      !print *, "Printing B_iaQ_bse_local of shape", SHAPE(B_iaQ_bse_local)
      !do i_print=1,SIZE(B_iaQ_bse_local,2)
      !   print *, B_iaQ_bse_local(:,:,1)
      !end do
      ! print *, "Printing B_abQ_bse_local of shape", SHAPE(B_abQ_bse_local)
      ! do k_print=1,SIZE(B_abQ_bse_local,3)
      !    do j_print=1,SIZE(B_abQ_bse_local,1)
      !       do i_print=1,SIZE(B_abQ_bse_local,2)
      !          print *,"(",i_print,j_print,k_print,")",B_abQ_bse_local(i_print,j_print,k_print)
      !          print *,"(",j_print,i_print,k_print,")",B_abQ_bse_local(j_print,i_print,k_print)
      !       end do
      !    end do
      ! end do
      !print *, "Printing B_bar_ijQ_bse_local of shape", SHAPE(B_bar_ijQ_bse_local)
      !do k_print=1,SIZE(B_bar_ijQ_bse_local,3)
      !    do j_print=1,SIZE(B_bar_ijQ_bse_local,1)
      !      do i_print=1,SIZE(B_bar_ijQ_bse_local,2)
      !         print *,"(",i_print,j_print,k_print,")",B_bar_ijQ_bse_local(i_print,j_print,k_print)
      !         print *,"(",j_print,i_print,k_print,")",B_bar_ijQ_bse_local(j_print,i_print,k_print)
      !     end do
      !   end do
      !end do

      i_iter = 1
      !DO i_iter = 1, num_davidson_iter

         print *, "Allocating Z_vec,AZ,BZ with dimensions (homo,virt,num_Z)",homo,virtual,num_Z_vectors
         !print *, 'you really enter here for i_iter', i_iter
         ALLOCATE (Z_vectors(homo, virtual, num_Z_vectors))
         Z_vectors = 0.0_dp

         ALLOCATE (AZ(homo, virtual, num_Z_vectors))
         AZ = 0.0_dp

         ALLOCATE(Z_vectors_reshaped(homo*virtual,num_Z_vectors))
         Z_vectors_reshaped = 0.0_dp

         ALLOCATE(AZ_reshaped(homo*virtual,num_Z_vectors))
         AZ_reshaped = 0.0_dp

         ALLOCATE (BZ(homo, virtual, num_Z_vectors))
         BZ = 0.0_dp

         local_RI_size = SIZE(B_iaQ_bse_local, 3)

         ALLOCATE (M_ia_tmp(homo, virtual))
         M_ia_tmp = 0.0_dp

         ALLOCATE (M_ji_tmp(homo, homo))
         M_ji_tmp = 0.0_dp

         ALLOCATE (RI_vector(local_RI_size, num_Z_vectors))
         RI_vector = 0.0_dp

         ALLOCATE(Subspace_eigenval(num_new_t))   
         Subspace_eigenval = 0.0_dp

         ALLOCATE(Subspace_eigenvec(num_Z_vectors,num_new_t))
         Subspace_eigenvec = 0.0_dp

         ALLOCATE(subspace_residuals_reshaped(homo*virtual,num_new_t))
         subspace_residuals_reshaped = 0.0_dp

         ALLOCATE(Subspace_add_dir(homo,virtual,num_new_t))
         Subspace_add_dir = 0.0_dp

         ALLOCATE(Preconditioner(homo,virtual,num_new_t))
         Preconditioner = 0.0_dp

         ALLOCATE(W_vectors(homo,virtual,num_Z_vectors+num_new_t))
         W_vectors = 0.0_dp



         if (i_iter == 1) then
            CALL initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)
         end if
         

         CALL compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, &
                         M_ia_tmp, RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, &
                         para_env)

         !CALL compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
         !                M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, &
         !                para_env)

         print *, 'finished AZ/BZ computation, AZ=',AZ


         !MG to check: Reshaping correct?
         AZ_reshaped = reshape(AZ,(/homo*virtual,num_Z_vectors/))
         Z_vectors_reshaped = reshape(Z_vectors,(/homo*virtual,num_Z_vectors/))

         ! Diagonalize M and extract smallest eigenvalues/corresponding eigenvector(Number: 1)
         CALL compute_diagonalize_ZAZ(AZ_reshaped, Z_vectors_reshaped,num_Z_vectors,homo,virtual,Subspace_eigenval, &
                                       Subspace_eigenvec,num_new_t)
         
         print *,"Eigenval",Subspace_eigenval
         print *,"Eigenvec",Subspace_eigenvec

         ! Compute residuals 
         CALL compute_residuals(AZ_reshaped, Z_vectors_reshaped,Subspace_eigenval,Subspace_eigenvec,Subspace_residuals_reshaped,&
                                 homo, virtual, num_new_t,num_Z_vectors)
         
         print *,"Residuals are"
         DO k_print=1,homo*virtual
            print *,Subspace_residuals_reshaped(k_print,:)
         END DO

         ! Calculate and add next orthonormal vector and update num_Z_vectors
         CALL compute_t(Preconditioner,homo,virtual,Subspace_residuals_reshaped,Subspace_eigenval,Eigenval, &
                        num_new_t,Subspace_add_dir)

         print *,"Subspace_add_dir"
         print *, Subspace_add_dir
         
         !Check dimensions and orthonormalize vector system
         if (num_Z_vectors + num_new_t .LE. num_max_z_space) then
            W_vectors(:,:,:num_Z_vectors) = Z_vectors(:,:,:)
            W_vectors(:,:,num_Z_vectors+num_new_t:) = Subspace_add_dir
            print *,"Added subspace_dir"
         else 
            print *,"****************************"
            print *,"Warning: Not yet implemented!"
            print *,"****************************"
            W_vectors = Z_vectors
         end if
         print *,"Z_vectors"
         DO k_print=1,num_Z_vectors
            print *,k_print
            print *,Z_vectors(:,:,k_print)
         END DO
         print *,"W_vectors"
         DO k_print=1,num_Z_vectors+num_new_t
            print *,k_print
            print *,W_vectors(:,:,k_print)
         END DO

         !Orthonorm:

         DEALLOCATE (AZ, BZ, Z_vectors, M_ia_tmp, M_ji_tmp, RI_vector, &
         Subspace_eigenval, Subspace_eigenvec,Subspace_residuals_reshaped,Subspace_add_dir)

         num_Z_vectors = num_Z_vectors + num_new_t
      !END DO
      

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE compute_t(Preconditioner,homo,virtual,Subspace_residuals_reshaped,Subspace_eigenval,Eigenval, &
                        num_new_t,Subspace_add_dir)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :,:)      :: Preconditioner,Subspace_add_dir
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_residuals_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_eigenval
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo,virtual,num_new_t

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_t'

      INTEGER                                            :: handle, mu_subspace, i_occ, a_virt
      REAL(KIND=dp)                                      :: prec_scalar
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_add_dir_reshaped

      ALLOCATE(Subspace_add_dir_reshaped(homo*virtual,num_new_t))

      CALL timeset(routineN,handle)
      print *,Subspace_residuals_reshaped
      DO mu_subspace=1,num_new_t
         DO i_occ = 1,homo
            Do a_virt = 1,virtual
               !MG to check: Indexorder and range of indices
                  prec_scalar = 1/(Subspace_eigenval(mu_subspace) - (Eigenval(a_virt + homo) - Eigenval(i_occ)))
                  print *, i_occ, a_virt, prec_scalar
                  Subspace_add_dir_reshaped((i_occ-1)*virtual+a_virt,mu_subspace) = prec_scalar * Subspace_residuals_reshaped((i_occ-1)*virtual+a_virt,mu_subspace)
            END DO
         END DO
      END DO
      print *,Subspace_add_dir_reshaped

      Subspace_add_dir = reshape(Subspace_add_dir_reshaped,(/homo,virtual,num_new_t/))

      DEALLOCATE(Subspace_add_dir_reshaped)
      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE compute_residuals(AZ_reshaped, Z_vectors_reshaped,Subspace_eigenval,Subspace_eigenvec,Subspace_residuals_reshaped, &
                                 homo, virtual, num_new_t,num_Z_vectors)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Subspace_eigenvec,AZ_reshaped, Z_vectors_reshaped,Subspace_residuals_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Subspace_eigenval
      INTEGER                                            :: homo,virtual,num_new_t,num_Z_vectors

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_residuals'

      INTEGER                                            :: handle,mu_subspace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)         :: Subspace_res_ev, Subspace_res_A

      CALL timeset(routineN,handle)

      ALLOCATE(Subspace_res_ev(homo*virtual,num_new_t))
      Subspace_res_ev = 0.0_dp

      ALLOCATE(Subspace_res_A(homo*virtual,num_new_t))
      Subspace_res_A = 0.0_dp

      !Compute all residuals in one loop, iterating over number of new/added t per iteration
      DO mu_subspace = 1,num_new_t
         
         CALL DGEMM("N", "N", homo*virtual, 1, num_Z_vectors, Subspace_eigenval(mu_subspace), Z_vectors_reshaped, homo*virtual, &
                       Subspace_eigenvec(:, mu_subspace), num_Z_vectors, 0.0_dp, Subspace_res_ev(:,mu_subspace), homo*virtual)

         CALL DGEMM("N", "N", homo*virtual, 1, num_Z_vectors, Subspace_eigenval(mu_subspace), AZ_reshaped, homo*virtual, &
                       Subspace_eigenvec(:, mu_subspace), num_Z_vectors, 0.0_dp, Subspace_res_A(:,mu_subspace), homo*virtual)

         Subspace_residuals_reshaped(:,mu_subspace) = Subspace_res_ev(:,mu_subspace) - Subspace_res_A(:,mu_subspace)

      END DO

      DEALLOCATE(Subspace_res_ev,Subspace_res_A)

      CALL timestop(handle)

   END SUBROUTINE


   SUBROUTINE compute_diagonalize_ZAZ(AZ_reshaped, Z_vectors_reshaped,num_Z_vectors,homo,virtual,Subspace_eigenval, &
                                       Subspace_eigenvec,num_new_t)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)                 :: AZ_reshaped, Z_vectors_reshaped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)                 :: Subspace_eigenvec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)                    :: Subspace_eigenval
      INTEGER,intent(in)                                          :: num_Z_vectors, homo, virtual,num_new_t

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_diagonalize_ZAZ'

      INTEGER                                               :: handle,  i_Z_vector, j_Z_vector, i_print, &
                                                               ZAZ_diag_info, LWORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)           :: ZAZ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)              :: WORK, eigenval_full_iter
      REAL(KIND=dp)                                         :: ZAZ_el

      ! interface lapack
      !    module subroutine DGEEV (jobvl,jobvr,N,A,LDA,WR,WI,VR,LDVR,WORK,LWORK,INFO)
      !       character,intent(in)       :: jobvl,jobvr
      !       integer,intent(in)         :: N,LDA,LDVR,LWORK
      !       integer, intent(out)       :: INFO
      !       REAL(kind=dp),intent(inout):: A(:,:)
      !       REAL(kind=dp),intent(out)  :: WR(:),WI(:),VR(:,:),WORK(:)
      !    end subroutine DGEEV
      ! end interface lapack

      CALL timeset(routineN,handle)
      
      ALLOCATE (ZAZ(num_Z_vectors,num_Z_vectors))
      ZAZ(:,:) = 0.0_dp
      
      ALLOCATE(eigenval_full_iter(num_Z_vectors))
      eigenval_full_iter = 0.0_dp


      !Flatten AZ and Z matrices of a certain j_Z_vector w.r.t. occ and virt indices
      !Multiply for each j_Z_vec and write into matrix of dim (num_Z_vec, num_Z_vec)
      DO i_Z_vector = 1, num_Z_vectors
         DO j_Z_vector = 1, num_Z_vectors
            ZAZ(j_Z_vector,i_Z_vector) = dot_product(Z_vectors_reshaped(:, j_Z_vector),AZ_reshaped(:, i_Z_vector))
         END DO
      END DO
      
      print *,"ZAZ matrix after loops:"
      DO i_print = 1, num_Z_vectors
         print *, ZAZ(i_print,:)
      END DO
      
      !Use hermitian/symmetric subroutine?
      LWORK = num_Z_vectors**2
      allocate(WORK(LWORK))
      !MG to check: Usage of symmetric routine okay? (Correct LWORK?)
      CALL DSYEV("V","U",num_Z_vectors,ZAZ,num_Z_vectors,eigenval_full_iter,WORK,LWORK,ZAZ_diag_info)

      if (ZAZ_diag_info /= 0) Then
         print *, "ZAZ could not be diagonalized successfully - Error:", ZAZ_diag_info
      else 
         print *,"Eigenvalues are", eigenval_full_iter
         print *,"Eigenvectors are", ZAZ
         print *,"WORK(1)=",WORK(1)
         print *,"ZAZ_diag_info",ZAZ_diag_info
      end if

      print *,"Shape of Subspace_eigenval",SHAPE(Subspace_eigenval)
      print *,"Shape of eigenval_full_iter",SHAPE(eigenval_full_iter)
      print *, "num_new_t",num_new_t
      print *, eigenval_full_iter(1:num_new_t)
      print *, Subspace_eigenval
      !print *, Subspace_eigenval(:)
      !print *, Subspace_eigenval(:num_new_t)
      Subspace_eigenval = eigenval_full_iter(1:num_new_t)
      print *,"Between eigenvalues and eigenvectors"
      Subspace_eigenvec(:,1:num_new_t) = ZAZ(:,1:num_new_t)
      print *,"After eigenvalues and eigenvectors"
      DEALLOCATE(WORK)
      DEALLOCATE(eigenval_full_iter)
      DEALLOCATE(ZAZ)

      CALL timestop(handle)
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param BZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param M_ji_tmp ...
!> \param RI_vector ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
                         M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ji_tmp, RI_vector
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i_Z_vector, LLL

      BZ(:, :, :) = 0.0_dp

      CALL compute_v_ia_jb_part(BZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                num_Z_vectors, homo, virtual)

      DO i_Z_vector = 1, num_Z_vectors

         DO LLL = 1, local_RI_size

            ! M_ji^P = sum_b Z_jb*B_bi^P
            CALL DGEMM("N", "T", homo, homo, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_iaQ_bse_local(:, :, LLL), homo, 0.0_dp, M_ji_tmp, homo)
            ! (BZ)_ia = sum_jP M_ij^P*B^bar_ja^P
            CALL DGEMM("T", "N", homo, virtual, homo, 1.0_dp, M_ji_tmp, homo, &
                       B_bar_iaQ_bse_local, homo, 1.0_dp, BZ(:, :, i_Z_vector), homo)

         END DO

      END DO

      ! we make the mp_sum to sum over all RI basis functions
      CALL mp_sum(BZ, para_env%group)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param M_ia_tmp ...
!> \param RI_vector ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, M_ia_tmp, &
                         RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_abQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ia_tmp, RI_vector
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_AZ'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, LLL
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      AZ(:, :, :) = 0.0_dp

      CALL compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                num_Z_vectors, homo, virtual)

      DO i_Z_vector = 1, num_Z_vectors

         ! JW TO DO: OMP PARALLELIZATION
         DO LLL = 1, local_RI_size

            ! M_ja^P = sum_b Z_jb*B_ba^P
            CALL DGEMM("N", "N", homo, virtual, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_abQ_bse_local(:, :, LLL), virtual, 0.0_dp, M_ia_tmp, homo)

            ! (AZ)_ia = sum_jP B_bar_ij^P*M_ja^P
            CALL DGEMM("N", "N", homo, virtual, homo, -1.0_dp, B_bar_ijQ_bse_local(:, :, LLL), homo, &
                       M_ia_tmp, homo, 1.0_dp, AZ(:, :, i_Z_vector), homo)

         END DO

      END DO

      ! we make the mp_sum to sum over all RI basis functions
      CALL mp_sum(AZ, para_env%group)

      ! add (e_a-e_i)*Z_ia
      DO i_occ = 1, homo
         DO a_virt = 1, virtual

            eigen_diff = Eigenval(a_virt + homo) - Eigenval(i_occ)

            AZ(i_occ, a_virt, :) = AZ(i_occ, a_virt, :) + Z_vectors(i_occ, a_virt, :)*eigen_diff

         END DO
      END DO


      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param RI_vector ...
!> \param local_RI_size ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                   num_Z_vectors, homo, virtual)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: AZ, Z_vectors, B_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: RI_vector
      INTEGER, INTENT(IN)                                :: local_RI_size, num_Z_vectors, homo, &
                                                            virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_v_ia_jb_part'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, LLL

      CALL timeset(routineN, handle)

      RI_vector = 0.0_dp

      ! v_P = sum_jb B_jb^P Z_jb
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual

                  RI_vector(LLL, i_Z_vector) = RI_vector(LLL, i_Z_vector) + &
                                               Z_vectors(i_occ, a_virt, i_Z_vector)* &
                                               B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO

      ! AZ = sum_P B_ia^P*v_P + ...
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual

                  AZ(i_occ, a_virt, i_Z_vector) = AZ(i_occ, a_virt, i_Z_vector) + &
                                                  RI_vector(LLL, i_Z_vector)* &
                                                  B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...Eigenval
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Z_vectors
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_Z_vectors, homo, virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initial_guess_Z_vectors'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, &
                                                            min_loc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigen_diff_ia

      CALL timeset(routineN, handle)

      ALLOCATE (eigen_diff_ia(homo, virtual))

      DO i_occ = 1, homo
         DO a_virt = 1, virtual
            eigen_diff_ia(i_occ, a_virt) = Eigenval(a_virt + homo) - Eigenval(i_occ)
         END DO
      END DO

      DO i_Z_vector = 1, num_Z_vectors

         min_loc = MINLOC(eigen_diff_ia)

         Z_vectors(min_loc(1), min_loc(2), i_Z_vector) = 1.0_dp

         eigen_diff_ia(min_loc(1), min_loc(2)) = 1.0E20_dp

      END DO

      DEALLOCATE (eigen_diff_ia)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
!> \param gd_array ...
!> \param color_sub ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE mult_B_with_W_and_fill_local_3c_arrays(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S, fm_mat_Q_static_bse, &
                                                     fm_mat_Q_static_bse_gemm, &
                                                     B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                     B_iaQ_bse_local, dimen_RI, homo, virtual, dimen_ia, &
                                                     gd_array, color_sub, para_env)

      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                                            fm_mat_S, fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual, dimen_ia
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mult_B_with_W_and_fill_local_3c_arrays'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_bar_ia_bse, &
                                                            fm_mat_S_bar_ij_bse, fm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_S_bar_ia_bse, fm_mat_S%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S, fm_mat_S_bar_ia_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ia_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_S_bar_ij_bse, fm_mat_S_ij_bse%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S_ij_bse, fm_mat_S_bar_ij_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ij_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_work, fm_mat_Q_static_bse_gemm%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_Q_static_bse_gemm, fm_mat_work)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      ! get info of fm_mat_Q_static_bse and compute ((1+Q(0))^-1-1)
      CALL cp_fm_get_info(matrix=fm_mat_Q_static_bse_gemm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) = fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) + 1.0_dp
            END IF
         END DO
      END DO

      ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q_static_bse_gemm, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      ! calculate [1+Q(i0)]^-1
      CALL cp_fm_cholesky_invert(fm_mat_Q_static_bse_gemm)
      ! symmetrize the result
      CALL cp_fm_upper_to_full(fm_mat_Q_static_bse_gemm, fm_mat_work)

      CALL cp_gemm(transa="N", transb="N", m=homo**2, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S_ij_bse, matrix_b=fm_mat_Q_static_bse, beta=0.0_dp, &
                   matrix_c=fm_mat_S_bar_ij_bse)

      ! fm_mat_S_bar_ia_bse has a different blacs_env as fm_mat_S_ij_bse since we take
      ! fm_mat_S from RPA. Therefore, we also need a different fm_mat_Q_static_bse_gemm
      CALL cp_gemm(transa="N", transb="N", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S, matrix_b=fm_mat_Q_static_bse_gemm, beta=0.0_dp, &
                   matrix_c=fm_mat_S_bar_ia_bse)

      CALL allocate_and_fill_local_array(B_iaQ_bse_local, fm_mat_S, gd_array, color_sub, homo, virtual, dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_iaQ_bse_local, fm_mat_S_bar_ia_bse, gd_array, color_sub, homo, virtual, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_ijQ_bse_local, fm_mat_S_bar_ij_bse, gd_array, color_sub, homo, homo, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_abQ_bse_local, fm_mat_S_ab_bse, gd_array, color_sub, virtual, virtual, &
                                         dimen_RI, para_env)

      CALL cp_fm_release(fm_mat_S_bar_ia_bse)
      CALL cp_fm_release(fm_mat_S_bar_ij_bse)

      CALL cp_fm_release(fm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param B_local ...
!> \param fm_mat_S ...
!> \param gd_array ...
!> \param color_sub ...
!> \param small_size ...
!> \param big_size ...
!> \param dimen_RI ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_local_array(B_local, fm_mat_S, gd_array, &
                                            color_sub, small_size, big_size, dimen_RI, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: B_local
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_mat_S
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub, small_size, big_size, dimen_RI
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_local_array'

      INTEGER :: combi_index, end_RI, handle, handle1, i_comm, i_entry, iiB, imepos, jjB, &
         level_big_size, level_small_size, ncol_local, nrow_local, num_comm_cycles, RI_index, &
         size_RI, start_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, mepos_from_RI_index, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp)                                      :: matrix_el
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      ALLOCATE (mepos_from_RI_index(dimen_RI))
      mepos_from_RI_index = 0

      DO imepos = 0, para_env%num_pe - 1

         CALL get_group_dist(gd_array, pos=imepos, starts=start_RI, ends=end_RI)

         mepos_from_RI_index(start_RI:end_RI) = imepos

      END DO

      ! color_sub is automatically the number of the process since every subgroup has only one MPI rank
      CALL get_group_dist(gd_array, color_sub, start_RI, end_RI, size_RI)

      ALLOCATE (B_local(small_size, big_size, 1:size_RI))

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))

      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      num_comm_cycles = 10

      ! communicate not all due to huge memory overhead, since for every number in fm_mat_S, we store
      ! three additional ones (RI index, first MO index, second MO index!!)
      DO i_comm = 0, num_comm_cycles - 1

         num_entries_send = 0
         num_entries_rec = 0

         ! loop over RI index to get the number of sent entries
         DO jjB = 1, ncol_local

            RI_index = col_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            num_entries_send(imepos) = num_entries_send(imepos) + nrow_local

         END DO

         CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

         ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
         ALLOCATE (buffer_send(0:para_env%num_pe - 1))

         ! allocate data message and corresponding indices
         DO imepos = 0, para_env%num_pe - 1

            ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
            buffer_rec(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
            buffer_send(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
            buffer_rec(imepos)%indx = 0

            ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
            buffer_send(imepos)%indx = 0

         END DO

         entry_counter(:) = 0

         ! loop over RI index for filling the send-buffer
         DO jjB = 1, ncol_local

            RI_index = col_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            DO iiB = 1, nrow_local

               combi_index = row_indices(iiB)
               level_small_size = MAX(1, combi_index - 1)/MAX(big_size,2) + 1
               level_big_size = combi_index - (level_small_size - 1)*big_size

               entry_counter(imepos) = entry_counter(imepos) + 1

               buffer_send(imepos)%msg(entry_counter(imepos)) = fm_mat_S%local_data(iiB, jjB)

               buffer_send(imepos)%indx(entry_counter(imepos), 1) = RI_index
               buffer_send(imepos)%indx(entry_counter(imepos), 2) = level_small_size
               buffer_send(imepos)%indx(entry_counter(imepos), 3) = level_big_size

            END DO

         END DO

         CALL timeset("BSE_comm_data", handle1)

         CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

         CALL timestop(handle1)

         ! fill B_local
         DO imepos = 0, para_env%num_pe - 1

            DO i_entry = 1, num_entries_rec(imepos)

               RI_index = buffer_rec(imepos)%indx(i_entry, 1) - start_RI + 1
               level_small_size = buffer_rec(imepos)%indx(i_entry, 2)
               level_big_size = buffer_rec(imepos)%indx(i_entry, 3)

               matrix_el = buffer_rec(imepos)%msg(i_entry)
               !MG to del
               print *,"level_small_size",level_small_size
               print *,"level_big_size",level_big_size
               print *,"RI_index",RI_index

               B_local(level_small_size, level_big_size, RI_index) = matrix_el

            END DO

         END DO

         DO imepos = 0, para_env%num_pe - 1
            DEALLOCATE (buffer_send(imepos)%msg)
            DEALLOCATE (buffer_send(imepos)%indx)
            DEALLOCATE (buffer_rec(imepos)%msg)
            DEALLOCATE (buffer_rec(imepos)%indx)
         END DO

         DEALLOCATE (buffer_rec, buffer_send)

      END DO

      DEALLOCATE (num_entries_send, num_entries_rec)

      DEALLOCATE (mepos_from_RI_index)

      DEALLOCATE (entry_counter, req_array)

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE deallocate_matrices_bse(fm_mat_S_ij_bse,fm_mat_S_ab_bse, &
                                      fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)

     TYPE(cp_fm_type), INTENT(INOUT), POINTER           :: fm_mat_S_ij_bse,fm_mat_S_ab_bse, &
                                                           fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm


     CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_bse'

     INTEGER                                            :: handle

     CALL timeset(routineN,handle)

     CALL cp_fm_release(fm_mat_S_ij_bse)
     CALL cp_fm_release(fm_mat_S_ab_bse)
     CALL cp_fm_release(fm_mat_Q_static_bse)
     CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

     CALL timestop(handle)
   END SUBROUTINE deallocate_matrices_bse

END MODULE bse
