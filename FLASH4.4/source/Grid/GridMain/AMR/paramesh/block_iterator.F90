!!****ih* source/Grid/GridMain/AMR/paramesh/block_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

module block_iterator

#include "Flash.h"
    use tree, ONLY : lnblocks, lrefine, lrefine_max
    use gr_specificData, ONLY : gr_oneBlock
    use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
    use Grid_interface, ONLY : Grid_getBlkIndexLimits
    use Grid_interface, ONLY : Grid_getBlkCornerID
    use Grid_interface, ONLY : Grid_blockMatch

    implicit none

#define CONTIGUOUS_POINTER pointer
#include "constants.h"
    private

    public :: destroy_iterator

    integer,parameter :: ndims=N_DIM
    integer,parameter :: amrex_real=kind(1.0)

    !!****ic* block_iterator/block_iterator_t
    !!
    !! NAME
    !!  block_iterator_t
    !!
    !!****
    type, public :: block_iterator_t
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
 
        this%nodetype = nodetype
        this%lev = INVALID_LEVEL 
        if (present(level)) then
            this%lev = level
        end if

        call this%first()
    end function init_iterator

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
        type(block_iterator_t), intent(INOUT) :: this

        call this%first()
    end subroutine destroy_iterator

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

        ! Search for the first valid block
        this%cur = 0
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
    logical function is_valid(this)
        use tree, ONLY : lnblocks

        class(block_iterator_t), intent(IN) :: this

        is_valid = (this%cur <= lnblocks)
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

        integer :: j = 0
  
        ! DEVNOTE: Move this check inside of Grid_blockMatch
        if (this%lev == INVALID_LEVEL) then
            do j = this%cur + 1, lnblocks
                if (Grid_blockMatch(j, this%nodetype)) EXIT
            end do
        else
            do j = this%cur + 1, lnblocks
                if (Grid_blockMatch(j, this%nodetype, this%lev)) EXIT
            end do
        end if

        this%cur = j
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
    subroutine blkMetaData(this, mData)
        use block_metadata, ONLY : block_metadata_t

        class(block_iterator_t), intent(IN)  :: this
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
    
end module block_iterator
