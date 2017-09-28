!!****ih* source/Grid/GridMain/AMR/Amrex/block_iterator
!!
!! This module is a facade pattern that maps the AMReX Fortran iterator onto 
!! the interface required presently by FLASH.
!!
!! Ideally, we will be able to use the AMReX iterator directly in the code 
!! and client code will gain access to it through implementation-specific 
!! code like Grid_getBlkIterator.
!!
!!****

module block_iterator

    use amrex_octree_module, ONLY : amrex_octree_iter, &
                                    amrex_octree_iter_build, &
                                    amrex_octree_iter_destroy

    implicit none

    private

#include "constants.h"
#include "Flash.h"

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
        type(amrex_octree_iter) :: oti
        integer                 :: level    = INVALID_LEVEL
        logical                 :: is_itor_valid = .FALSE.
    contains
        procedure, public :: is_valid
        procedure, public :: first
        procedure, public :: next
        procedure, public :: blkMetaData
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
        final             :: destroy_iterator
#endif
    end type block_iterator_t

    interface block_iterator_t
        procedure :: init_iterator
    end interface block_iterator_t

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
    function init_iterator(nodetype, level) result(this)
        integer, intent(IN)           :: nodetype
        integer, intent(IN), optional :: level
        type(block_iterator_t)        :: this

        if (present(level)) then
            this%level = level
        end if

        ! DEVNOTE: the AMReX iterator is not built based on nodetype.
        ! It appears that we get leaves every time.

        ! Initial iterator is not primed.  Advance to first compatible block.
        call amrex_octree_iter_build(this%oti)
        call this%next()
    end function init_iterator

#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
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
    subroutine destroy_iterator(this)
        type(block_iterator_t), intent(INOUT) :: this

        call amrex_octree_iter_destroy(this%oti)
    end subroutine destroy_iterator
#endif

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

        ! reset to before first valid block
        call this%oti%clear()   !method added to amrex_octree_iter - KW 2017-08-23
        this%is_itor_valid = .FALSE.  !DEV: ??
        ! Initial iterator is not primed.  Advance to first compatible block.
        call this%next()
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

        ans = this%is_itor_valid
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

        this%is_itor_valid = this%oti%next()

        if (this%level /= INVALID_LEVEL) then
            ! Search for leaves on given level
            do while (this%is_itor_valid)
                ! oti has 0-based level indexing, while this has 1-based
                if (this%oti%level() == (this%level - 1)) then
                    exit
                end if

                this%is_itor_valid = this%oti%next()
            end do
        end if
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
        use amrex_box_module,     ONLY : amrex_box

        use block_metadata,       ONLY : block_metadata_t
        use gr_physicalMultifabs, ONLY : unk

        class(block_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: blockDesc

        integer         :: n_guards(MDIM) = 0
        type(amrex_box) :: box
   
        box = this%oti%box()

        ! Block descriptor provides FLASH-compliant 1-based level index set,
        ! but AMReX uses 0-based index set.
        blockDesc%level             = this%oti%level() + 1
        blockDesc%grid_index        = this%oti%grid_index()
        ! Block descriptor provides FLASH-compliant 1-based cell index set,
        ! but AMReX uses 0-based index set.
        blockDesc%limits(LOW,  :) = 1
        blockDesc%limits(HIGH, :) = 1
        blockDesc%limits(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
        blockDesc%limits(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1

        ! DEVNOTE: KW says that box with GC available through newer AMReX
        ! fortran interface.
        n_guards(:) = 0
        ! Multifab arrays are 0-based (AMReX) instead of 1-based(FLASH)
        n_guards(1:NDIM) = unk(blockDesc%level-1)%nghost() 
        blockDesc%limitsGC(LOW,  :) = 1 
        blockDesc%limitsGC(HIGH, :) = 1
        blockDesc%limitsGC(LOW,  1:NDIM) =   blockDesc%limits(LOW,  1:NDIM) &
                                           - n_guards(1:NDIM)
        blockDesc%limitsGC(HIGH, 1:NDIM) =   blockDesc%limits(HIGH, 1:NDIM) &
                                           + n_guards(1:NDIM)

        blockDesc%localLimits(LOW,  :)      = 1
        blockDesc%localLimits(HIGH, :)      = 1
        blockDesc%localLimits(LOW,  1:NDIM) = NGUARD + 1
        blockDesc%localLimits(HIGH, 1:NDIM) =   blockDesc%limits(HIGH, 1:NDIM) &
                                              - blockDesc%limits(LOW,  1:NDIM) + NGUARD + 1
        blockDesc%localLimitsGC(LOW,  :) = 1
        blockDesc%localLimitsGC(HIGH, :) = 1
        blockDesc%localLimitsGC(HIGH, 1:NDIM) =   blockDesc%limitsGC(HIGH,1:NDIM) &
                                                - blockDesc%limitsGC(LOW, 1:NDIM) &
                                                + 1

    end subroutine blkMetaData
 
end module block_iterator

