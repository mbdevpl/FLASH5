!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection2/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(       real :: initData,
!!                       integer(IN) :: tileDesc) 
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

subroutine Simulation_initBlock(initData, tileDesc)
    use flash_tile, ONLY : flash_tile_t 

    implicit none

    real,                          pointer :: initData(:, :, :, :)
    type(flash_tile_t), intent(IN)         :: tileDesc

    integer :: i, j, k, var

    associate(lo => tileDesc%limits(LOW,  :), &
              hi => tileDesc%limits(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        initData(i, j, k, var) = 1.0
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

