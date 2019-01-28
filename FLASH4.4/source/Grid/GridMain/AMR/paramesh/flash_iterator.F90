!!****ih* source/Grid/GridMain/AMR/paramesh/flash_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"
#include "Flash.h"
#include "constants.h"

module flash_iterator
    use tree, ONLY : lnblocks

    implicit none

    private
    public :: build_iterator, destroy_iterator

    !!****ic* flash_iterator/flash_iterator_t
    !!
    !! NAME
    !!  flash_iterator
    !!
    !!****
    type, public :: flash_iterator_t
        integer :: cur      = 0
        integer :: nodetype = LEAF 
        integer :: lev      = INVALID_LEVEL
    contains
        procedure, public :: isValid
        procedure, public :: next
        procedure, public :: currentTile
    end type flash_iterator_t

contains

    !!****im* flash_iterator/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(flash_iterator_t(OUT)    :: itor,
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
    !!             Refer to the documentation for the Paramesh version of
    !!             gr_getBlkIterator for more details.
    !!  level    - iterate only over all blocks/tiles of the correct nodetype
    !!             that are located at this level of refinement.  Refer to the
    !!             documentation for the Paramesh version of gr_getBlkIterator
    !!             for more details.  A level value of UNSPEC_LEVEL is equivalent
    !!             to omitting this optional argument.
    !!  tiling   - an optional optimization hint.  Tiling is not implemented for
    !!             Paramesh and therefore this hint is ignored.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, nodetype, level, tiling, tileSize)
        type(flash_iterator_t), intent(OUT) :: itor
        integer,                intent(IN)  :: nodetype
        integer,                intent(IN)  :: level
        logical,                intent(IN)  :: tiling
        integer,                intent(IN)  :: tileSize(1:MDIM)

        itor%nodetype = nodetype
        itor%lev = level

        itor%cur = 0
        call itor%next()
    end subroutine build_iterator

    !!****im* flash_iterator_t/destroy_iterator
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
        type(flash_iterator_t), intent(INOUT) :: itor

        itor%cur = lnblocks + 1
    end subroutine destroy_iterator

    !!****m* flash_iterator_t/isValid
    !!
    !! NAME
    !!  isValid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%isValid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block
    !!
    !!****
    logical function isValid(this)
        class(flash_iterator_t), intent(IN) :: this

        isValid = (this%cur <= lnblocks)
    end function isValid

    !!****m* flash_iterator_t/next
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

        class(flash_iterator_t), intent(INOUT) :: this

        integer :: j

        if (this%lev == UNSPEC_LEVEL) then
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

    !!****m* flash_iterator_t/blkMetaData
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
    subroutine currentTile(this, tileDesc)
        use gr_specificData,            ONLY : gr_oneBlock
        use flash_tile,                 ONLY : flash_tile_t 
        use tree,                       ONLY : lrefine, &
                                               lrefine_max

        class(flash_iterator_t), intent(IN)  :: this
        type(flash_tile_t),      intent(OUT) :: tileDesc

        integer :: cornerID(1:MDIM)
        integer :: blkLim(LOW:HIGH, 1:MDIM)

        tileDesc%id = this%cur 
        tileDesc%level = lrefine(tileDesc%id)

        tileDesc%cid = gr_oneBlock(tileDesc%id)%cornerID
        tileDesc%stride = 2**(lrefine_max - tileDesc%level)

        associate(lo      => tileDesc%limits(LOW, :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGrown => tileDesc%grownLimits(LOW, :), &
                  hiGrown => tileDesc%grownLimits(HIGH, :), &
                  loBlkGC => tileDesc%blkLimitsGC(LOW, :), &
                  hiBlkGC => tileDesc%blkLimitsGC(HIGH, :), &
                  blkId   => tileDesc%id, &
                  cid     => tileDesc%cid)
            cornerID = (cid - 1) / 2**(lrefine_max-lrefine(blkID)) + 1

            ! Get limits with block-local indexing of cells
            blkLim(LOW,  IAXIS) = GRID_ILO
            blkLim(HIGH, IAXIS) = GRID_IHI
            blkLim(LOW,  JAXIS) = GRID_JLO
            blkLim(HIGH, JAXIS) = GRID_JHI
            blkLim(LOW,  KAXIS) = GRID_KLO
            blkLim(HIGH, KAXIS) = GRID_KHI

            ! Convert to global indexing of cells
            lo(:)   = blkLim(LOW,  :) - 1 + cornerID(:)
            hi(:)   = blkLim(HIGH, :) - 1 + cornerID(:)
            lo(1:NDIM)   = lo(1:NDIM) - NGUARD
            hi(1:NDIM)   = hi(1:NDIM) - NGUARD

            loBlkGC(:) = lo(:)
            hiBlkGC(:) = hi(:)
            loBlkGC(1:NDIM) = loBlkGC(1:NDIM) - NGUARD
            hiBlkGC(1:NDIM) = hiBlkGC(1:NDIM) + NGUARD

            ! Since there is no tiling with Paramesh, the grown tile
            ! is just the block + GC halo
            loGrown(:) = loBlkGC(:)
            hiGrown(:) = hiBlkGC(:)
        end associate
    end subroutine currentTile

end module flash_iterator

