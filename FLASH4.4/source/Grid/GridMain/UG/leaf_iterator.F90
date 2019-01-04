!!****ih* source/Grid/GridMain/AMR/paramesh/leaf_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

module leaf_iterator

#include "Flash.h"

    implicit none

#define CONTIGUOUS_POINTER pointer
#include "constants.h"
    private

    public :: build_iterator, destroy_iterator

    integer,parameter :: ndims=N_DIM
    integer,parameter :: amrex_real=kind(1.0)

    !!****ic* leaf_iterator/leaf_iterator_t
    !!
    !! NAME
    !!  leaf_iterator_t
    !!
    !!****
    type, public :: leaf_iterator_t
        integer :: cur         = 1
        integer :: lev         = INVALID_LEVEL
        integer :: cellIdxBase = -2
        ! different variants for cell index numbering:
        !  1 = normal FLASH convention:     per block, leftmost guard cell = 1
        !  0 = zero-based FLASH convention: per block, leftmost guard cell = 0
        ! -1 = global convention for a refinement level: leftmost guard cell = 1
        ! -2 = global convention for a refinement level: leftmost inner cell = 1
    contains
        procedure, public :: first
        procedure, public :: is_valid
        procedure, public :: next
        procedure, public :: blkMetaData
    end type leaf_iterator_t

contains

    !!****im* leaf_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(leaf_iterator_t(OUT)  :: itor,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of leaf blocks or
    !!  tiles within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!  
    !! ARGUMENTS
    !!  itor     - the iterator
    !!  level    - iterate only over leaf blocks/tiles located at this level of
    !!             refinement.
    !!             A level value of UNSPEC_LEVEL is equivalent to omitting
    !!             this optional argument.
    !!  tiling   - an optional optimization hint.  Tiling is not implemented for
    !!             Paramesh and therefore this hint is ignored.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, level, tiling)
        type(leaf_iterator_t), intent(OUT)          :: itor
        integer,               intent(IN), optional :: level
        logical,               intent(IN), optional :: tiling

        itor%lev = 1

        call itor%first()
    end subroutine build_iterator

    !!****im* leaf_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Destroy given iterator.
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
        type(leaf_iterator_t), intent(INOUT) :: itor

        call itor%first()
    end subroutine destroy_iterator

    !!****m* leaf_iterator_t/first
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
        class(leaf_iterator_t), intent(INOUT) :: this

        ! Search for the first valid block
        this%cur = 0
        call this%next()
    end subroutine first
 
    !!****m* leaf_iterator_t/is_valid
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
    logical function is_valid(this)

        class(leaf_iterator_t), intent(IN) :: this

        is_valid = (this%cur <= 1)
    end function is_valid

    !!****m* leaf_iterator_t/next
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
        class(leaf_iterator_t), intent(INOUT) :: this

        this%cur = 2
    end subroutine next

    !!****m* leaf_iterator_t/blkMetaData
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
    subroutine blkMetaData(this, mData)
        use block_metadata,             ONLY : block_metadata_t
        use Grid_getBlkIndexLimits_mod, ONLY : Grid_getBlkIndexLimits
        use Grid_data, ONLY : 

        class(leaf_iterator_t), intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: mData

        integer, dimension(MDIM)           :: cornerID
        integer, dimension(LOW:HIGH, MDIM) :: blkLim, blkLimGC
        
        mData%id = this%cur 
        mData%level = 1

        ! DEVNOTE: For this to work, all metadata must be determined here
        ! without using any of the Grid_get* subroutines that take
        ! blockId as an argument
        mData%cid = gr_blkCornerID
        mData%stride = 1

        associate(lo    => mData%limits(LOW, :), &
                  hi    => mData%limits(HIGH, :), &
                  loGC  => mData%limitsGC(LOW, :), &
                  hiGC  => mData%limitsGC(HIGH, :), &
                  blkId => mData%id, &
                  cid   => mData%cid)
            call Grid_getBlkIndexLimits(1, blkLim, blkLimGC)
            mdata%localLimits   = blkLim
            mdata%localLimitsGC = blkLimGC
            lo(:) = blkLim(LOW, :)
            hi(:) = blkLim(HIGH, :)
            loGC(:) = blkLimGC(LOW, :)
            hiGC(:) = blkLimGC(HIGH, :)
        end associate
    end subroutine blkMetaData
    
end module leaf_iterator
