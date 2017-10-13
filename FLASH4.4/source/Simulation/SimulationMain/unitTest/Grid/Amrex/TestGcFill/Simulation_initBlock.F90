!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestInit/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a composit number which is a combination
!!  of the block number and the indices of the cell
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

#include "constants.h"
#include "Flash.h"
 
subroutine Simulation_initBlock(initData, block)
    use Grid_interface, ONLY : Grid_getSingleCellCoords
    use block_metadata, ONLY : block_metadata_t, bmd_print
    
    implicit none
    
    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: block

    integer :: i, j, k, var

    integer :: idx(1:MDIM)
    real    :: coords(1:MDIM)

    ! Set data only in interiors so that GC are set to zero
    initData(:, :, :, :) = 0.0d0
    associate(lo => block%limits(LOW,  :), &
              hi => block%limits(HIGH, :))
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    idx = [i - lo(IAXIS) + 1, &
                           j - lo(JAXIS) + 1, &
                           k - lo(KAXIS) + 1]
                    call Grid_getSingleCellCoords(idx, block, &
                                                  CENTER, INTERIOR, coords)
                    do var=UNK_VARS_BEGIN, UNK_VARS_END
                        initData(i, j, k, var) = &
                             DBLE((coords(IAXIS) + coords(JAXIS)) * var)
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

