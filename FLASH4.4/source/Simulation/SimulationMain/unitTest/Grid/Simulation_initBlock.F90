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

subroutine Simulation_initBlock(solnData,blockDesc)
  use Grid_interface, ONLY :  Grid_putPointData, Grid_getBlkPtr, Grid_getCellCoords, Grid_getSingleCellVol
  use Grid_interface, ONLY :  Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshME
  use block_metadata, ONLY : block_metadata_t

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  real,dimension(:,:,:,:),pointer :: solnData,sfacex,sfacey
  type(block_metadata_t), intent(in) :: blockDesc

  integer :: i, j, k, var,blockID
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: grd,size,point
  real,dimension(:),allocatable::xLeft,xRight,xCenter,yCenter,yLeft,yRight

  logical,parameter :: gcell=.true.
  real :: areaFactor_1, areaFactor_2, areaFactor_3


  blkLimits = blockDesc%limits
  blkLimitsGC = blockDesc%limitsGC
  blockID=blockDesc%id
  grd(:)=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
  size=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
#define TEST_USE_FACEVARS
#ifdef TEST_USE_FACEVARS
  call Grid_getBlkPtr(blockDesc,sfacex,FACEX)
  call Grid_getBlkPtr(blockDesc,sfacey,FACEY)
#endif
  allocate(xLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xRight(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xCenter(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))

  allocate(yLeft(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yRight(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yCenter(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))

  call Grid_getCellCoords(IAXIS,blockDesc,LEFT_EDGE,gcell,xLeft,size(IAXIS))
  call Grid_getCellCoords(IAXIS,blockDesc,RIGHT_EDGE,gcell,xRight,size(IAXIS))
  call Grid_getCellCoords(IAXIS,blockDesc,LEFT_EDGE,gcell,xCenter,size(IAXIS))

  call Grid_getCellCoords(JAXIS,blockDesc,LEFT_EDGE,gcell,yLeft,size(JAXIS))
  call Grid_getCellCoords(JAXIS,blockDesc,RIGHT_EDGE,gcell,yRight,size(JAXIS))
  call Grid_getCellCoords(JAXIS,blockDesc,LEFT_EDGE,gcell,yCenter,size(JAXIS))
  
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           point(IAXIS)=i; point(JAXIS)=j; point(KAXIS)=k
            areaFactor_1 = xLeft(i)
            areaFactor_2 = yRight(j) - yLeft(j)
            areaFactor_3 = 2.*PI
#ifdef TEST_USE_FACEVARS
            sfacex(1,i,j,k) = abs(areaFactor_1 * areaFactor_2 * areaFactor_3)
#endif
            call Grid_getSingleCellVol(blockDesc, point, solnData(1,i,j,k))
        end do
#ifdef TEST_USE_FACEVARS
        associate (i => blkLimitsGC(HIGH,IAXIS))
          sfacex(1,i+1,j,k)= abs(xRight(i) * areaFactor_2 * areaFactor_3)
        end associate
#endif
     end do
  end do
!!$  blk=(gr_MeshMe*100.0+blockID)*100000000.0
!!$  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
!!$     kind = k*1000000.0
!!$    do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
!!$       jind = j*10000.0
!!$       do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
!!$          iind = i*100.0
!!$          do var=UNK_VARS_BEGIN,UNK_VARS_END
!!$             val = blk+kind+jind+iind+var
!!$             solnData(var,i,j,k)=val
!!$             sfacex(1,i,j,k)=val
!!$             sfacey(1,i,j,k)=val
!!$          end do
!!$          sfacex(1,blkLimitsGC(HIGH,IAXIS)+1,j,k)=blk+kind+jind+iind+100+1
!!$       enddo
!!$    enddo
!!$ enddo
!!$
!!$ j=blkLimitsGC(HIGH,JAXIS)+1
!!$ jind = j*10000.0
!!$ do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
!!$    kind = k*1000000.0
!!$    do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
!!$       iind = i*100.0
!!$       do var=UNK_VARS_BEGIN,UNK_VARS_END
!!$          val = blk+kind+jind+iind+var
!!$          sfacey(1,i,j,k)=blk+kind+jind+iind+100+1
!!$       enddo
!!$    end do
!!$ end do

#ifdef TEST_USE_FACEVARS
  call Grid_ReleaseBlkPtr(blockDesc,sfacex,FACEX)
  call Grid_ReleaseBlkPtr(blockDesc,sfacey,FACEY)
#endif
 return
end subroutine Simulation_initBlock
