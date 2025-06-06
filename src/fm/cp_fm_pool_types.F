!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief pool for for elements that are retained and released
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_fm_pool_types
   USE cp_fm_struct, ONLY: cp_fm_struct_release, &
                           cp_fm_struct_retain, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: cp_fm_create, &
                          cp_fm_p_type, &
                          cp_fm_release, &
                          cp_fm_type
   USE cp_linked_list_fm, ONLY: cp_sll_fm_dealloc, &
                                cp_sll_fm_get_first_el, &
                                cp_sll_fm_insert_el, &
                                cp_sll_fm_next, &
                                cp_sll_fm_rm_first_el, &
                                cp_sll_fm_type
   USE cp_log_handling, ONLY: cp_to_string
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_pool_types'

   PUBLIC :: cp_fm_pool_type, cp_fm_pool_p_type
   PUBLIC :: fm_pool_create, fm_pool_retain, &
             fm_pool_release, &
             fm_pool_create_fm, fm_pool_give_back_fm, &
             fm_pool_get_el_struct
   PUBLIC :: fm_pools_dealloc, &
             fm_pools_create_fm_vect, &
             fm_pools_give_back_fm_vect
!***

! **************************************************************************************************
!> \brief represent a pool of elements with the same structure
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param el_struct the structure of the elements stored in this pool
!> \param cache linked list with the elements in the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE cp_fm_pool_type
      PRIVATE
      INTEGER :: ref_count = -1
      TYPE(cp_fm_struct_type), POINTER :: el_struct => NULL()
      TYPE(cp_sll_fm_type), POINTER :: cache => NULL()
   END TYPE cp_fm_pool_type

! **************************************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE cp_fm_pool_p_type
      TYPE(cp_fm_pool_type), POINTER :: pool => NULL()
   END TYPE cp_fm_pool_p_type

   INTERFACE fm_pools_create_fm_vect
      MODULE PROCEDURE fm_pools_create_fm_m1_p_type_pointer
      MODULE PROCEDURE fm_pools_create_fm_m1_p_type_alloc
      MODULE PROCEDURE fm_pools_create_fm_m1_array_pointer
      MODULE PROCEDURE fm_pools_create_fm_m1_array_alloc
   END INTERFACE

   INTERFACE fm_pools_give_back_fm_vect
      MODULE PROCEDURE fm_pools_give_back_fm_m1_p_type_pointer
      MODULE PROCEDURE fm_pools_give_back_fm_m1_p_type_alloc
      MODULE PROCEDURE fm_pools_give_back_fm_m1_array_pointer
      MODULE PROCEDURE fm_pools_give_back_fm_m1_array_alloc
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief creates a pool of elements
!> \param pool the pool to create
!> \param el_struct the structure of the elements that are stored in
!>        this pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pool_create(pool, el_struct)
      TYPE(cp_fm_pool_type), POINTER                     :: pool
      TYPE(cp_fm_struct_type), TARGET                    :: el_struct

      ALLOCATE (pool)
      pool%el_struct => el_struct
      CALL cp_fm_struct_retain(pool%el_struct)
      pool%ref_count = 1

   END SUBROUTINE fm_pool_create

! **************************************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pool_retain(pool)
      TYPE(cp_fm_pool_type), INTENT(INOUT)               :: pool

      CPASSERT(pool%ref_count > 0)

      pool%ref_count = pool%ref_count + 1
   END SUBROUTINE fm_pool_retain

! **************************************************************************************************
!> \brief deallocates all the cached elements
!> \param pool the pool to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pool_flush_cache(pool)
      TYPE(cp_fm_pool_type), INTENT(IN)                  :: pool

      TYPE(cp_fm_type), POINTER                          :: el_att
      TYPE(cp_sll_fm_type), POINTER                      :: iterator

      iterator => pool%cache
      DO
         IF (.NOT. cp_sll_fm_next(iterator, el_att=el_att)) EXIT
         CALL cp_fm_release(el_att)
         DEALLOCATE (el_att)
         NULLIFY (el_att)
      END DO
      CALL cp_sll_fm_dealloc(pool%cache)
   END SUBROUTINE fm_pool_flush_cache

! **************************************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pool_release(pool)
      TYPE(cp_fm_pool_type), POINTER                     :: pool

      IF (ASSOCIATED(pool)) THEN
         CPASSERT(pool%ref_count > 0)
         pool%ref_count = pool%ref_count - 1
         IF (pool%ref_count == 0) THEN
            pool%ref_count = 1
            CALL fm_pool_flush_cache(pool)
            CALL cp_fm_struct_release(pool%el_struct)
            pool%ref_count = 0

            DEALLOCATE (pool)
         END IF
      END IF
      NULLIFY (pool)
   END SUBROUTINE fm_pool_release

! **************************************************************************************************
!> \brief returns an element, allocating it if none is in the pool
!> \param pool the pool from where you get the element
!> \param element will contain the new element
!>\param name the name for the new matrix (optional)
!> \param name ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pool_create_fm(pool, element, &
                                name)
      TYPE(cp_fm_pool_type), INTENT(IN)                  :: pool
      TYPE(cp_fm_type), INTENT(OUT)                      :: element
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: name

      TYPE(cp_fm_type), POINTER                          :: el

      NULLIFY (el)
      IF (ASSOCIATED(pool%cache)) THEN
         el => cp_sll_fm_get_first_el(pool%cache)
         CALL cp_sll_fm_rm_first_el(pool%cache)
      END IF
      IF (ASSOCIATED(el)) THEN
         element = el
         DEALLOCATE (el)
      ELSE
         CALL cp_fm_create(element, matrix_struct=pool%el_struct)
      END IF

      IF (PRESENT(name)) THEN
         element%name = name
      ELSE
         element%name = "tmp"
      END IF

   END SUBROUTINE fm_pool_create_fm

! **************************************************************************************************
!> \brief returns the element to the pool
!> \param pool the pool where to cache the element
!> \param element the element to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      transfers the ownership of the element to the pool
!>      (it is as if you had called cp_fm_release)
!>      Accept give_backs of non associated elements?
! **************************************************************************************************
   SUBROUTINE fm_pool_give_back_fm(pool, element)
      TYPE(cp_fm_pool_type), INTENT(IN)                  :: pool
      TYPE(cp_fm_type), INTENT(INOUT)                    :: element

      IF (.NOT. ASSOCIATED(pool%el_struct, element%matrix_struct)) THEN
         CALL cp_fm_release(element)
      ELSE
         BLOCK
            TYPE(cp_fm_type), POINTER :: el
            ALLOCATE (el)
            el = element
            CALL cp_sll_fm_insert_el(pool%cache, el=el)
            NULLIFY (element%matrix_struct, element%local_data, element%local_data_sp)
         END BLOCK
      END IF
   END SUBROUTINE fm_pool_give_back_fm

! **************************************************************************************************
!> \brief returns the structure of the elements in this pool
!> \param pool the pool you are interested in
!> \return ...
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   FUNCTION fm_pool_get_el_struct(pool) RESULT(res)
      TYPE(cp_fm_pool_type), INTENT(IN)                  :: pool
      TYPE(cp_fm_struct_type), POINTER                   :: res

      res => pool%el_struct
   END FUNCTION fm_pool_get_el_struct

!================== pools ================

! **************************************************************************************************
!> \brief shallow copy of an array of pools (retains each pool)
!> \param source_pools the pools to copy
!> \param target_pools will contains the new pools
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pools_copy(source_pools, target_pools)
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: source_pools, target_pools

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(source_pools))
      ALLOCATE (target_pools(SIZE(source_pools)))
      DO i = 1, SIZE(source_pools)
         target_pools(i)%pool => source_pools(i)%pool
         CALL fm_pool_retain(source_pools(i)%pool)
      END DO
   END SUBROUTINE fm_pools_copy

! **************************************************************************************************
!> \brief deallocate an array of pools (releasing each pool)
!> \param pools the pools to release
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fm_pools_dealloc(pools)
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: pools

      INTEGER                                            :: i

      IF (ASSOCIATED(pools)) THEN
         DO i = 1, SIZE(pools)
            CALL fm_pool_release(pools(i)%pool)
         END DO
         DEALLOCATE (pools)
      END IF
   END SUBROUTINE fm_pools_dealloc

   #:mute
      #:set types = [("cp_fm_type", "array", ""), ("cp_fm_p_type", "p_type","%matrix")]
      #:set attributes = [("ALLOCATABLE", "alloc", "ALLOCATED"), ("POINTER", "pointer", "ASSOCIATED")]
   #:endmute

   #:for typename, shortname, appendix in types
      #:for attr, shortattr, create in attributes
! **************************************************************************************************
!> \brief Returns a vector with an element from each pool
!> \param pools the pools to create the elements from
!> \param elements will contain the vector of elements
!> \param name the name for the new matrixes (optional)
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE fm_pools_create_fm_m1_${shortname}$_${shortattr}$ (pools, elements, &
                                                                       name)
            TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(IN)  :: pools
            TYPE(${typename}$), DIMENSION(:), ${attr}$          :: elements
            CHARACTER(len=*), INTENT(in), OPTIONAL             :: name

            INTEGER                                            :: i
            TYPE(cp_fm_pool_type), POINTER                     :: pool

            NULLIFY (pool)

            ALLOCATE (elements(SIZE(pools)))
            DO i = 1, SIZE(pools)
               pool => pools(i)%pool
               #:if typename=="cp_fm_p_type"
                  ALLOCATE (elements(i)%matrix)
               #:endif
               IF (PRESENT(name)) THEN
                  CALL fm_pool_create_fm(pool, elements(i) ${appendix}$, &
                                         name=name//"-"//ADJUSTL(cp_to_string(i)))
               ELSE
                  CALL fm_pool_create_fm(pool, elements(i) ${appendix}$)
               END IF

            END DO

         END SUBROUTINE fm_pools_create_fm_m1_${shortname}$_${shortattr}$

! **************************************************************************************************
!> \brief returns a vector to the pools. The vector is deallocated
!>      (like cp_fm_vect_dealloc)
!> \param pools the pool where to give back the vector
!> \param elements the vector of elements to give back
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      accept unassociated vect?
! **************************************************************************************************
         SUBROUTINE fm_pools_give_back_fm_m1_${shortname}$_${shortattr}$ (pools, elements)
            TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(IN)  :: pools
            TYPE(${typename}$), DIMENSION(:), ${attr}$         :: elements

            INTEGER                                            :: i

            IF (${create}$ (elements)) THEN
               CPASSERT(SIZE(pools) == SIZE(elements))
               DO i = 1, SIZE(pools)
                  CALL fm_pool_give_back_fm(pools(i)%pool, &
                                            elements(i) ${appendix}$)
                  #:if typename == "cp_fm_p_type"
                     DEALLOCATE (elements(i)%matrix)
                  #:endif
               END DO
               DEALLOCATE (elements)
               #:if attr == "POINTER"
                  NULLIFY (elements)
               #:endif
            END IF
         END SUBROUTINE fm_pools_give_back_fm_m1_${shortname}$_${shortattr}$
      #:endfor
   #:endfor

END MODULE cp_fm_pool_types
