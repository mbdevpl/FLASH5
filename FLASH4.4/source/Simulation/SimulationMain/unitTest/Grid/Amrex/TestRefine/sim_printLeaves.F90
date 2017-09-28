subroutine sim_printLeaves(title, block_count)
    use gr_amrexInterface, ONLY : gr_getFinestLevel
    use block_iterator,    ONLY : block_iterator_t
    use block_metadata,    ONLY : block_metadata_t

    implicit none

#include "constants.h"

    character(*), intent(IN)    :: title
    integer,      intent(INOUT) :: block_count(4)

    type(block_iterator_t) :: itor
    type(block_metadata_t) :: blockDesc

    integer :: lev
    integer :: level
    integer :: finest_level

    block_count(:) = 0

    write(*,'(A)') title
    write(*,'(A)') "-----------------------------------"
    call gr_getFinestLevel(finest_level)
    do level = 1, finest_level 
        write(*,'(A,I2)') "Leaf blocks at level ", level
        itor = block_iterator_t(LEAF, level=level)
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)
            write(*,'(A,I3,A,I3,A,I3,A,I3,A)') "     From (", &
                     blockDesc%limits(LOW, 1), ", ", &
                     blockDesc%limits(LOW, 2), ") to (", &
                     blockDesc%limits(HIGH, 1), ", ", &
                     blockDesc%limits(HIGH, 2), ")"
    
            block_count(level) = block_count(level) + 1

            call itor%next()
        end do
    end do
end subroutine sim_printLeaves

