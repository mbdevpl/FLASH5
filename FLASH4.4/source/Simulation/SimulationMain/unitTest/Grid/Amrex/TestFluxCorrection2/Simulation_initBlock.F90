!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection2/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(       real :: initData,
!!                       integer(IN) :: block) 
!!
!! DESCRIPTION
!!  Initializes the cell-centered data to one.  For the purposes of this test, 
!!  the initial value of the cell-centered data is unimportant.
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

    associate(lo => block%limitsGC(LOW,  :), &
              hi => block%limitsGC(HIGH, :))
        do               k = lo(KAXIS), hi(KAXIS)
            do           j = lo(JAXIS), hi(JAXIS)
                do       i = lo(IAXIS), hi(IAXIS)
                    do var = UNK_VARS_BEGIN, UNK_VARS_END
                        initData(i, j, k, var) = 1.0
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

