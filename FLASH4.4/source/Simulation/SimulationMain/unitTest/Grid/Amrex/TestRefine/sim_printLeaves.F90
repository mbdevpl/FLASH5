subroutine sim_printLeaves(title)
    use gr_amrexInterface, ONLY : gr_getFinestLevel
    use block_iterator,    ONLY : block_iterator_t
    use block_metadata,    ONLY : block_metadata_t
    use Simulation_data,   ONLY : leaves, &
                                  MIN_REFINE_LEVEL

    implicit none

#include "constants.h"

    character(*), intent(IN)    :: title

    type(block_iterator_t) :: itor
    type(block_metadata_t) :: blockDesc

    integer :: level, j
    integer :: finest_level

    write(*,'(A)') title
    write(*,'(A)') "-----------------------------------"
    call gr_getFinestLevel(finest_level)
    do level = MIN_REFINE_LEVEL, finest_level 
        if (.NOT. allocated(leaves(level)%blocks))   CYCLE
        
        associate(blks => leaves(level)%blocks)
            write(*,'(A,I2)') "Leaf blocks at level ", level
            do j = 1, SIZE(leaves(level)%blocks, 1)
                write(*,'(A,I3,A,I3,A,I3,A,I3,A)') "     From (", &
                         blks(j, 1), ", ", &
                         blks(j, 2), ") to (", &
                         blks(j, 3), ", ", &
                         blks(j, 4), ")"
            end do
        end associate

    end do
end subroutine sim_printLeaves

