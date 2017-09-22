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
    use block_metadata, ONLY : block_metadata_t
    use sim_interface,  ONLY : sim_writeDataPoints

    implicit none

    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: block

#include "constants.h"
#include "Flash.h"

    real    :: points(2, 2)
    real    :: values(2)

    initData(:, :, :, :) = 0.0d0

    points(:, :) = 0.0d0
    points(1, :) = [0.16, 0.67]
    points(2, :) = [0.11, 0.38]
    values(:) = 0.0d0
    values(1) = 3.0d0
    values(2) = 1.0d0

    call sim_writeDataPoints(initData, block, points, values)
end subroutine Simulation_initBlock

