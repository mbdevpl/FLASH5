!!****ih* source/Grid/GridMain/AMR/paramesh/gr_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

module gr_iterator

#include "Flash.h"
    use tree, ONLY : lnblocks, lrefine, lrefine_max
    use gr_specificData, ONLY : gr_oneBlock

    implicit none

#define CONTIGUOUS_POINTER pointer
#include "constants.h"
    private

    public :: build_iterator, destroy_iterator

    integer,parameter :: ndims=N_DIM
    integer,parameter :: amrex_real=kind(1.0)

    !!****ic* gr_iterator/gr_iterator_t
    !!
    !! NAME
    !!  gr_iterator_t
    !!
    !!****
    type, public :: gr_iterator_t
        integer :: cur         = 1
        integer :: nodetype    = LEAF 
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
    end type gr_iterator_t

contains

    !!****im* gr_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(gr_iterator_t(OUT)    :: itor,
    !!                 integer(IN)           :: nodetype,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!  
    !! ARGUMENTS
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).
    !!             Refer to the documetation for the Paramesh version of 
    !!             gr_getBlkIterator for more details.
    !!  level    - iterate only over all blocks/tiles of the correct nodetype
    !!             that are located at this level of refinement.  Refer to the
    !!             documetation for the Paramesh version of gr_getBlkIterator 
    !!             for more details.
    !!  tiling   - an optional optimization hint.  Tiling is not implemented for
    !!             Paramesh and therefore this hint is ignored.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, nodetype, level, tiling)
        type(gr_iterator_t),    intent(OUT)          :: itor
        integer,                intent(IN)           :: nodetype
        integer,                intent(IN), optional :: level
        logical,                intent(IN), optional :: tiling

        itor%nodetype = nodetype
        itor%lev = INVALID_LEVEL 
        if (present(level)) then
            itor%lev = level
        end if

        call itor%first()
    end subroutine build_iterator

    !!****im* gr_iterator_t/destroy_iterator
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
        type(gr_iterator_t), intent(INOUT) :: itor

        call itor%first()
    end subroutine destroy_iterator

    !!****m* gr_iterator_t/first
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
        class(gr_iterator_t), intent(INOUT) :: this

        ! Search for the first valid block
        this%cur = 0
        call this%next()
    end subroutine first
 
    !!****m* gr_iterator_t/is_valid
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
        use tree, ONLY : lnblocks

        class(gr_iterator_t), intent(IN) :: this

        is_valid = (this%cur <= lnblocks)
    end function is_valid

    !!****m* gr_iterator_t/next
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
        use gr_parameshInterface, ONLY : gr_blockMatch

        class(gr_iterator_t), intent(INOUT) :: this

        integer :: j = 0
  
        if (this%lev == INVALID_LEVEL) then
            ! No level given at creation
            do j = this%cur + 1, lnblocks
                if (gr_blockMatch(j, this%nodetype)) EXIT
            end do
        else
            do j = this%cur + 1, lnblocks
                if (gr_blockMatch(j, this%nodetype, this%lev)) EXIT
            end do
        end if

        this%cur = j
    end subroutine next

    !!****m* gr_iterator_t/blkMetaData
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

        class(gr_iterator_t),    intent(IN)  :: this
        type(block_metadata_t),  intent(OUT) :: mData

        integer, dimension(MDIM)           :: cornerID
        integer, dimension(LOW:HIGH, MDIM) :: blkLim, blkLimGC
        
        mData%id = this%cur 
        mData%level = lrefine(mData%id)

        ! DEVNOTE: For this to work, all metadata must be determined here
        ! without using any of the Grid_get* subroutines that take
        ! blockId as an argument
        mData%cid = gr_oneBlock(mData%id)%cornerID
        mData%stride = 2**(lrefine_max - mData%level)

        associate(lo    => mData%limits(LOW, :), &
                  hi    => mData%limits(HIGH, :), &
                  loGC  => mData%limitsGC(LOW, :), &
                  hiGC  => mData%limitsGC(HIGH, :), &
                  blkId => mData%id, &
                  cid   => mData%cid)
            call Grid_getBlkIndexLimits(blkID, blkLim, blkLimGC)
            mdata%localLimits   = blkLim
            mdata%localLimitsGC = blkLimGC
            lo(:) = blkLim(LOW, :)
            hi(:) = blkLim(HIGH, :)
            loGC(:) = blkLimGC(LOW, :)
            hiGC(:) = blkLimGC(HIGH, :)
            if (this%cellIdxBase == -1) then
               cornerID = (cid - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:)   = lo(:)   - 1 + cornerID(:)
               hi(:)   = hi(:)   - 1 + cornerID(:)
               loGC(:) = loGC(:) - 1 + cornerID(:)
               hiGC(:) = hiGC(:) - 1 + cornerID(:)
            else if (this%cellIdxBase == -2) then
               cornerID = (cid - 1) / 2**(lrefine_max-lrefine(blkID)) + 1
               lo(:)   = lo(:)   - 1 + cornerID(:)
               hi(:)   = hi(:)   - 1 + cornerID(:)
               loGC(:) = loGC(:) - 1 + cornerID(:)
               hiGC(:) = hiGC(:) - 1 + cornerID(:)
               lo(1:ndims)   = lo(1:ndims)   - NGUARD
               hi(1:ndims)   = hi(1:ndims)   - NGUARD
               loGC(1:ndims) = loGC(1:ndims) - NGUARD
               hiGC(1:ndims) = hiGC(1:ndims) - NGUARD
            else if (this%cellIdxBase == 0) then
               lo(:)   = lo(:)   - 1
               hi(:)   = hi(:)   - 1
               loGC(:) = loGC(:) - 1
               hiGC(:) = hiGC(:) - 1
            end if
        end associate
    end subroutine blkMetaData
    
end module gr_iterator
