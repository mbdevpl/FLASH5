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

    integer :: i, j, k
    integer :: var

    integer :: point(2)

    associate(lo => block%limits(LOW,  :), &
              hi => block%limits(HIGH, :))
        initData(:, :, :, :) = 0.0d0

        point = [3, 11] * (2**(block%level - 1))
        if (      (lo(1) <= point(1)) .AND. (point(1) <= hi(1)) &
            .AND. (lo(2) <= point(2)) .AND. (point(2) <= hi(2))) then
            initData(point(1), point(2), 1, 1) = 2.0d0
        end if
        
        point = [3, 6]
        if (      (lo(1) <= point(1)) .AND. (point(1) <= hi(1)) &
            .AND. (lo(2) <= point(2)) .AND. (point(2) <= hi(2))) then
            initData(point(1), point(2), 1, 1) = 1.0d0
        end if
    end associate
end subroutine Simulation_initBlock

