!!****ih* source/Grid/GridMain/UG/flash_iterator
!!
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "constants.h"
#include "FortranLangFeatures.fh"

module flash_iterator

#include "Flash.h"

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* flash_iterator/flash_iterator_t
    !!
    !! NAME
    !!  flash_iterator_t
    !!
    !!****
    type, public :: flash_iterator_t
        integer :: cur = 2
        integer :: lev = INVALID_LEVEL
    contains
        procedure, public :: isValid
        procedure, public :: next
        procedure, public :: currentTile 
    end type flash_iterator_t

contains

    !!****im* flash_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(flash_iterator_t(OUT)  :: itor,
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
    subroutine build_iterator(itor, nodetype, level, tiling)
        type(flash_iterator_t), intent(OUT)          :: itor
        integer,                intent(IN)           :: nodetype
        integer,                intent(IN), optional :: level
        logical,                intent(IN), optional :: tiling

        itor%cur = 1
        itor%lev = 1

        if (present(level)) then
            itor%lev = level
            if      (itor%lev == UNSPEC_LEVEL) then
                itor%lev = 1    
            else if (itor%lev /= 1) then
                itor%cur = 2
            end if
        end if

        if ((nodetype /= LEAF) .AND. (nodetype /= ALL_BLKS)) then
            call Driver_abortFlash("[flash_iterator]: Unsupported nodetype")
        end if
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

        itor%cur = 2
        itor%lev = INVALID_LEVEL
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

        isValid = (this%cur == 1)
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
        class(flash_iterator_t), intent(INOUT) :: this

        this%cur = 2
    end subroutine next

    !!****m* flash_iterator_t/currentTile
    !!
    !! NAME
    !!  currentTile 
    !!
    !! SYNPOSIS
    !!  call itor%currentTile(flash_tile_t(OUT) : block)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine currentTile(this, tileDesc)
        use flash_tile,       ONLY : flash_tile_t
        use Grid_data,        ONLY : gr_blkCornerID, &
                                     gr_ilo, gr_ihi, &
                                     gr_jlo, gr_jhi, &
                                     gr_klo, gr_khi, &
                                     gr_iloGc, gr_ihiGc, &
                                     gr_jloGc, gr_jhiGc, &
                                     gr_kloGc, gr_khiGc
        use Driver_interface, ONLY : Driver_abortFlash

        class(flash_iterator_t), intent(IN)  :: this
        type(flash_tile_t),      intent(OUT) :: tileDesc

        if (.NOT. this%isValid()) then
            call Driver_abortFlash("[currentTile] No current tile")
        end if

        tileDesc%id = this%cur 
        tileDesc%level = this%lev
        tileDesc%cid = gr_blkCornerID
        tileDesc%stride = 1

        associate(lo      => tileDesc%limits(LOW, :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGrown => tileDesc%grownLimits(LOW, :), &
                  hiGrown => tileDesc%grownLimits(HIGH, :), &
                  loBlkGC => tileDesc%blkLimitsGC(LOW, :), &
                  hiBlkGC => tileDesc%blkLimitsGC(HIGH, :))
            lo(1:MDIM)   = [gr_ilo,   gr_jlo,   gr_klo]
            hi(1:MDIM)   = [gr_ihi,   gr_jhi,   gr_khi]
            loBlkGC(1:MDIM) = [gr_iloGc, gr_jloGc, gr_kloGc]
            hiBlkGC(1:MDIM) = [gr_ihiGc, gr_jhiGc, gr_khiGc]

            ! Since there is no tiling with UG, the grown tile
            ! is just the block + GC halo
            loGrown(:) = loBlkGC(:)
            hiGrown(:) = hiBlkGC(:)
        end associate
    end subroutine currentTile

end module flash_iterator

