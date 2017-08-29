!!****ih* source/Grid/GridMain/AMR/Amrex/block_iterator
!!
!! This module is a facade pattern that maps the AMReX Fortran iterator onto 
!! the interface required presently by FLASH.
!!
!! Ideally, we will be able to use the AMReX iterator directly in the code 
!! and client code will gain access to it through implementation-specific 
!! code like Grid_getBlkIterator.
!!
!! This is a variant that uses type(block_1lev_iterator_t) as underlying iterator.
!!****

#include "FortranLangFeatures.fh"
#include "constants.h"

module block_iterator

    use amrex_multifab_module, ONLY : amrex_mfiter, &
                                    amrex_mfiter_build, &
                                    amrex_mfiter_destroy
    use block_1lev_iterator !, ONLY : block_1lev_iterator_t

    implicit none

    private

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
        type(block_1lev_iterator_t),allocatable :: li(:)
        integer                 :: first_level   = INVALID_LEVEL
        integer                 :: last_level    = INVALID_LEVEL
        integer                 :: level    = INVALID_LEVEL
        logical                 :: isValid = .FALSE.
        real,POINTER            :: fp(:,:,:,:)
    contains
        procedure, public :: is_valid
        procedure, public :: first
        procedure, public :: next
        procedure, public :: blkMetaData
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
        final             :: destroy_iterator
#else
        procedure         :: destroy_iterator
#endif
    end type block_iterator_t

    interface block_iterator_t
        procedure :: init_iterator
        procedure :: init_iterator_mfa
    end interface block_iterator_t

    interface block_iterator_destroy
       procedure :: destroy_iterator
    end interface block_iterator_destroy

contains

    !!****im* block_iterator_t/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !! SYNOPOSIS
    !!  block_iterator_t itor = block_iterator_t(integer(IN)         :: nodetype,
    !!                                           level(IN), optional :: level)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks
    !!  within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block.
    !!
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
    !!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
    !!             iterate only over blocks located at this level of 
    !!             octree structure
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    function init_iterator_mfa(nodetype, mfArray, level, tiling) result(this)
      use amrex_multifab_module, ONLY : amrex_multifab

        type(block_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        type(amrex_multifab),intent(IN),CONTIGUOUS :: mfArray(:)
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

        integer :: l, first, last
        logical :: v

        if (present(level)) then
            first = level
            last = level
         else
            first = 1
            last = size(mfArray)
        end if
 
        allocate( this%li (first : last) )

        this%first_level = first
        this%last_level = last
        this%level = first

!!$        print*,'block_iterator_build: about to build 1lev iterators for this=',this%isValid,this%level,allocated(this%li)

        do l=first,last
!!$           call amrex_mfiter_build(this%li(l),mfArray(l),tiling=tiling)
            this%li(l) = block_1lev_iterator_t(nodetype, mfArray(l),l,tiling=tiling)
!!$            call this%li( l )%first()
            v = this%li( l )%is_valid()
            if (v .AND. .NOT. this%isValid) then
               this%isValid = .TRUE.
               this%level   = l
            end if
        end do

        if (.NOT. this%isValid) then
           call this%destroy_iterator()
        end if

!!$        print*,'block_iterator_build: done building 1lev iterators for this=',this%isValid,this%level,allocated(this%li)
!!$        call this%first()
      end function init_iterator_mfa

    function init_iterator(nodetype, level, tiling) result(this)
      use amrex_multifab_module, ONLY : amrex_multifab
      use gr_amrextData

        type(block_iterator_t)        :: this
        integer, intent(IN)           :: nodetype
        integer, intent(IN), optional :: level
        logical, intent(IN), optional :: tiling

!!$        type(amrex_multifab),POINTER :: mfArray(:)
!!$        type(amrex_multifab),POINTER :: mf

        this = init_iterator_mfa(nodetype, gr_amrextUnkMFs, level, tiling)
    end function init_iterator

!#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
    !!****im* block_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Called automatically
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(this)
      class (block_iterator_t), intent(INOUT) :: this

        integer :: l

        if (allocated(this%li)) then
           do l = this%first_level, this%last_level

              call this%li( l )%destroy_iterator()

           end do
           deallocate(this%li)
        end if
        this%isValid = .FALSE.

    end subroutine destroy_iterator
!#endif

    !!****m* block_iterator_t/first
    !!
    !! NAME
    !!  first
    !!
    !! SYNPOSIS
    !!  call itor%first() 
    !!
    !! DESCRIPTION
    !!  Reset iterator to the initial block managed by process
    !!
    !!****
    subroutine first(this)
        class(block_iterator_t), intent(INOUT) :: this

        integer :: l
        logical :: v

        call Driver_abortFlash('block_iterator: Attempting first(), not implemented!')
        print*,'block_iterator%first: about to do 1lev%first''s on this=',this%isValid,this%level,allocated(this%li)

        do l = this%first_level, this%last_level
           call this%li( l )%first()
        end do

        print*,'block_iterator%first: done 1lev%first''s on this=',this%isValid,this%level,allocated(this%li)

        if (this%first_level .LE. this%last_level) then
           l = this%first_level

           v = this%li( l )%is_valid()

           do while (l .LE. this%last_level .AND. .NOT. v)
              l = l+1
              call this%li( l )%first()
              v = this%li( l )%is_valid()
           end do


           this%level = l
           this%isValid = v

        end if



    end subroutine first
 
    !!****m* block_iterator_t/is_valid
    !!
    !! NAME
    !!  is_valid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%is_valid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block
    !!
    !!****
    function is_valid(this) result(ans)
        class(block_iterator_t), intent(IN) :: this
        logical :: ans

        ans = this%isValid
    end function is_valid

    !!****m* block_iterator_t/next
    !!
    !! NAME
    !!  next
    !!
    !! SYNPOSIS
    !!  call itor%next()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next block managed by process and that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        class(block_iterator_t), intent(INOUT) :: this

        integer :: l
        logical :: v

!!$        print*,'block_iterator%next: about to do 1lev%next on this=',this%isValid,this%level,allocated(this%li)

        l = this%level

        call this%li( l )%next()
        v = this%li( l )%is_valid()

        do while (l .LT. this%last_level .AND. .NOT. v)
           l = l+1
!!$           call this%li( l )%first()   ! not necessary now!
           v = this%li( l )%is_valid()
        end do

        this%level = l
        this%isValid = v

    end subroutine next

    !!****m* block_iterator_t/blkMetaData
    !!
    !! NAME
    !!  blkMetaData 
    !!
    !! SYNPOSIS
    !!  call itor%blkMetaData(block_metadata_t(OUT) : block)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine blkMetaData(this, blockDesc)
        use amrex_box_module, ONLY : amrex_box
        use block_metadata, ONLY : block_metadata_t
        use tree,           ONLY : lrefine_max

        class(block_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: blockDesc

        type(amrex_box) :: box, fabbox
       
        box = this%li( this%level )%tilebox()
        fabbox=this%li(this%level )%fabbox()

        ! TODO: Determine if box contains GC or not and finalize limits/limitsGC
!!$        blockDesc%grid_index        = this%oti%grid_index()
        blockDesc%level             = this%level
        blockDesc%limits(LOW, :)    = box%lo
        blockDesc%limits(HIGH, :)   = box%hi
        blockDesc%limitsGC(LOW, :)  = fabbox%lo
        blockDesc%limitsGC(HIGH, :) = fabbox%hi

        blockDesc%localLimits(LOW, :)   = blockDesc%limits(LOW, :)   - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimits(HIGH, :)  = blockDesc%limits(HIGH, :)  - blockDesc%limitsGC(LOW, :) + 1
        blockDesc%localLimitsGC(LOW, :) = 1
        blockDesc%localLimitsGC(HIGH, :)= blockDesc%limitsGC(HIGH, :)- blockDesc%limitsGC(LOW, :) + 1

        blockDesc%fp => this%li( this%level )%dataPtr()

        blockDesc%id = -999
        blockDesc%stride = 2**(lrefine_max - blockDesc%level)
        blockDesc%cid = (box%lo - 1) * blockDesc%stride + 1

    end subroutine blkMetaData
 
end module block_iterator

