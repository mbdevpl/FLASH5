#include "constants.h"

subroutine sim_collectLeaves
    use Grid_interface,       ONLY : Grid_getBlkIterator, &
                                     Grid_releaseBlkIterator
    use gr_amrexInterface,    ONLY : gr_getFinestLevel
    use block_iterator,       ONLY : block_iterator_t 
    use block_metadata,       ONLY : block_metadata_t 
    use Simulation_data,      ONLY : leaves, &
                                     MIN_REFINE_LEVEL, MAX_REFINE_LEVEL

    type(block_iterator_t) :: itor
    type(block_metadata_t) :: blockDesc

    logical :: gridChanged
    integer :: finest_level
    integer :: block_count
    integer :: lev, j

    ! Regenerate leaf block data structure
    call gr_getFinestLevel(finest_level)
    do lev = MIN_REFINE_LEVEL, MAX_REFINE_LEVEL
        block_count = 0
        call Grid_getBlkIterator(itor, LEAF, level=lev)
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)
 
            block_count = block_count + 1

            call itor%next()
        end do
        call Grid_releaseBlkIterator(itor)

        if (allocated(leaves(lev)%blocks)) then
            deallocate(leaves(lev)%blocks)
        end if

        if (block_count > 0) then
            allocate(leaves(lev)%blocks(block_count, 4))
        end if
    end do

    ! Populate leaf block data structure
    do lev = MIN_REFINE_LEVEL, finest_level 
        call Grid_getBlkIterator(itor, LEAF, level=lev)

        j = 1
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)

            associate(lo => blockDesc%limits(LOW,  :), &
                      hi => blockDesc%limits(HIGH, :))
                leaves(lev)%blocks(j, :) = [lo(IAXIS), lo(JAXIS), &
                                            hi(IAXIS), hi(JAXIS)]
            end associate

            j = j + 1
            call itor%next()
        end do
        
        call Grid_releaseBlkIterator(itor)
    end do
end subroutine sim_collectLeaves

