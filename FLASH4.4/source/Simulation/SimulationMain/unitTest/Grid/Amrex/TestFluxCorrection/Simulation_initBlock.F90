!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(       real :: initData,
!!                       integer(IN) :: block) 
!!
!! DESCRIPTION
!!  Initializes the data to the 1-based level at which the block exists
!! 
!! ARGUMENTS
!!  initData - the cell-centered data structure to which data is written
!!  block - index of block whose cell-centered data is to be initialized 
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Simulation_initBlock(initData, block)
    use block_metadata, ONLY : block_metadata_t

    implicit none

    real,                   pointer    :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN) :: block

    integer :: i, j, k, var

    ! Initialize data.  The values are not important for this test
    associate(lo => block%limits(LOW,  :), &
              hi => block%limits(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        initData(i, j, k, var) = block%level
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

