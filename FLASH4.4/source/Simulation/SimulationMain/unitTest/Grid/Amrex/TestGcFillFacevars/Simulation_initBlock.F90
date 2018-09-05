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
    use Grid_interface, ONLY : Grid_getSingleCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr
    use block_metadata, ONLY : block_metadata_t, bmd_print
    
    implicit none
    
    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: block

    integer :: i, j, k, var

    integer :: idx(1:MDIM)
    real    :: coords(1:MDIM)
    real, dimension(:,:,:,:), pointer :: faceDataX, faceDataY, faceDataZ

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

!Set face centererd data in interiors
#if NFACE_VARS>0
        !Set faceDataX
        call Grid_getBlkPtr(block,faceDataX,FACEX)
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)+1
                    idx = [i - lo(IAXIS) + 1, &
                           j - lo(JAXIS) + 1, &
                           k - lo(KAXIS) + 1]
                    call Grid_getSingleCellCoords(idx, block, &
                                                  LEFT_EDGE, INTERIOR, coords)  !Coord of LEFT_EDGE
                    do var=1,NFACE_VARS
                        faceDataX(i, j, k, var) = &
                             DBLE((coords(IAXIS) + coords(JAXIS)) * var)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block,faceDataX,FACEX)
#if NDIM>1
        !Set faceDataY
        call Grid_getBlkPtr(block,faceDataY,FACEY)
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)+1
                do i = lo(IAXIS), hi(IAXIS)
                    idx = [i - lo(IAXIS) + 1, &
                           j - lo(JAXIS) + 1, &
                           k - lo(KAXIS) + 1]
                    call Grid_getSingleCellCoords(idx, block, &
                                                  LEFT_EDGE, INTERIOR, coords)  !Coord of LEFT_EDGE
                    do var=1,NFACE_VARS
                        faceDataY(i, j, k, var) = &
                             DBLE((coords(IAXIS) + coords(JAXIS)) * var)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block,faceDataY,FACEY)
#endif
#if NDIM>2
        !Set faceDataZ
        call Grid_getBlkPtr(block,faceDataZ,FACEZ)
        do         k = lo(KAXIS), hi(KAXIS)+1
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    idx = [i - lo(IAXIS) + 1, &
                           j - lo(JAXIS) + 1, &
                           k - lo(KAXIS) + 1]
                    call Grid_getSingleCellCoords(idx, block, &
                                                  LEFT_EDGE, INTERIOR, coords)  !Coord of LEFT_EDGE
                    do var=1,NFACE_VARS
                        faceDataZ(i, j, k, var) = &
                             DBLE((coords(IAXIS) + coords(JAXIS)) * var)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block,faceDataZ,FACEZ)
#endif

#endif
    end associate
end subroutine Simulation_initBlock

