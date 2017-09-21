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

subroutine Simulation_initBlock(initData, block)
    use block_metadata, ONLY : block_metadata_t, bmd_print
    
    implicit none
    
    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: block

#include "constants.h"
#include "Flash.h"
  
    integer :: i = 1
    integer :: j = 1
    integer :: k = 1
    integer :: var = 1

    associate(lo => block%limitsGC(LOW,  :), &
              hi => block%limitsGC(HIGH, :))
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    do var=UNK_VARS_BEGIN, UNK_VARS_END
                        initData(i, j, k, var) = 1.1d0 * var
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

