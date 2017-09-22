subroutine sim_printLeaves(title)
    use amrex_amrcore_module,  ONLY : amrex_get_finest_level
    
    use block_iterator,        ONLY : block_iterator_t
    use block_metadata,        ONLY : block_metadata_t

    implicit none

#include "constants.h"

    character(*), intent(IN) :: title

    type(block_iterator_t) :: itor
    type(block_metadata_t) :: blockDesc

    integer :: lev
    integer :: level
    integer :: finest_level
    
    write(*,'(A)') title
    write(*,'(A)') "-------------------------------------------------------"
    finest_level = amrex_get_finest_level() + 1
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

            call itor%next()
        end do
    end do
end subroutine sim_printLeaves

