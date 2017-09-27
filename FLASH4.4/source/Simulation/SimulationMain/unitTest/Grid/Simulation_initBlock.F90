!!****if* source/Simulation/SimulationMain/unitTest/Grid/Simulation_initBlock
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

subroutine Simulation_initBlock(solnData,block)
  use Grid_interface, ONLY :  Grid_putPointData
  use Grid_data, ONLY : gr_meshME
  use block_metadata, ONLY : block_metadata_t

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block

  integer :: i, j, k, var
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: grd,axis
  
  real :: blk,kind,jind,iind,val
  integer :: blockID

  blkLimits = block%limits
  blkLimitsGC = block%limitsGC
  blockID=block%id
  grd(:)=blkLimits(LOW,:)-blkLimitsGC(LOW,:)

  blk=(gr_MeshMe*100.0+blockID)*100000000.0
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     kind = k*1000000.0
     axis(KAXIS)=k
    do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
       jind = j*10000.0
       axis(JAXIS)=j
       do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
          iind = i*100.0
          axis(IAXIS)=i
          do var=UNK_VARS_BEGIN,UNK_VARS_END
             val = blk+kind+jind+iind+var
             solnData(var,i,j,k)=val
          end do
       enddo
    enddo
 enddo
  
 return
end subroutine Simulation_initBlock
