!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFaceVar/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!! DESCRIPTION
!!  Initializes the Grid with a composit number which is a combination
!!  of the block number and the indices of the cell
!! 
!! ARGUMENTS
!!  blockId -          the blockId to update
!!
!!***

#include "constants.h"
#include "Flash.h"
  
subroutine Simulation_initBlock(initData, block)
    use block_metadata, ONLY : block_metadata_t
    use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
    
    implicit none

    real,                   intent(IN), pointer :: initData(:, :, :, :)
    type(block_metadata_t), intent(IN)          :: block

    integer :: i, j, k, var

    real, pointer :: faceDataX(:,:,:,:) => null()
    real, pointer :: faceDataY(:,:,:,:) => null()
    real, pointer :: faceDataZ(:,:,:,:) => null()

#if NFACE_VARS>0
    ! Set face centererd data in interiors and on boundaries
    associate(lo => block%limits(LOW,  :), &
              hi => block%limits(HIGH, :))
        call Grid_getBlkPtr(block, faceDataX, FACEX)
        do           var = 1, NFACE_VARS
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)+1
                        faceDataX(i, j, k, var) = REAL(1.3*i*i)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block, faceDataX, FACEX)
#if NDIM>1
        call Grid_getBlkPtr(block, faceDataY, FACEY)
        do           var = 1, NFACE_VARS
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)+1
                    do i = lo(IAXIS), hi(IAXIS)
                        faceDataY(i, j, k, var) = REAL(i*i + 1.2*j*j)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block, faceDataY, FACEY)
#endif
#if NDIM>2
        call Grid_getBlkPtr(block,faceDataZ,FACEZ)
        do           var = 1, NFACE_VARS
            do         k = lo(KAXIS), hi(KAXIS)+1
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        faceDataZ(i, j, k, var) = REAL(i*i + j*j + 1.1*k*k)
                    end do
                end do
            end do
        end do
        call Grid_releaseBlkPtr(block,faceDataZ,FACEZ)
#endif
    end associate
#endif

end subroutine Simulation_initBlock

