#include "constants.h"

subroutine sim_collectLeaves
    use Grid_interface,       ONLY : Grid_getTileIterator, &
                                     Grid_releaseTileIterator
    use gr_amrexInterface,    ONLY : gr_getFinestLevel
    use flash_iterator,       ONLY : flash_iterator_t 
    use flash_tile,           ONLY : flash_tile_t
    use Simulation_data,      ONLY : leaves, &
                                     MIN_REFINE_LEVEL, MAX_REFINE_LEVEL

    type(flash_iterator_t) :: itor
    type(flash_tile_t)     :: blockDesc

    logical :: gridChanged
    integer :: finest_level
    integer :: block_count
    integer :: lev, j

    ! Regenerate leaf block data structure
    call gr_getFinestLevel(finest_level)
    do lev = MIN_REFINE_LEVEL, MAX_REFINE_LEVEL
        block_count = 0
        call Grid_getTileIterator(itor, LEAF, level=lev, tiling=.FALSE.)
        do while (itor%isValid())
            call itor%currentTile(blockDesc)
 
            block_count = block_count + 1

            call itor%next()
        end do
        call Grid_releaseTileIterator(itor)

        if (allocated(leaves(lev)%blocks)) then
            deallocate(leaves(lev)%blocks)
        end if

        if (block_count > 0) then
            allocate(leaves(lev)%blocks(block_count, 4))
        end if
    end do

    ! Populate leaf block data structure
    do lev = MIN_REFINE_LEVEL, finest_level 
        call Grid_getTileIterator(itor, LEAF, level=lev, tiling=.FALSE.)

        j = 1
        do while (itor%isValid())
            call itor%currentTile(blockDesc)

            associate(lo => blockDesc%limits(LOW,  :), &
                      hi => blockDesc%limits(HIGH, :))
                leaves(lev)%blocks(j, :) = [lo(IAXIS), lo(JAXIS), &
                                            hi(IAXIS), hi(JAXIS)]
            end associate

            j = j + 1
            call itor%next()
        end do

        call Grid_releaseTileIterator(itor)
    end do
end subroutine sim_collectLeaves

