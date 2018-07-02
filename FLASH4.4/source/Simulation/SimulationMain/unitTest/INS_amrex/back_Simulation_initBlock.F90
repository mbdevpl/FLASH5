!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! 
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData,block)

!  use Simulation_data
  use Simulation_data, ONLY :sim_xMin,sim_xMax,sim_yMin,sim_yMax,sim_zMin,sim_zMax
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block
  integer :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)
  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord
  integer :: sizeX,sizeY,sizeZ

  real :: Lx, Ly, Lz, xi, yi, zi, Phi_ijk, F_ijk

  real :: del(3)

  real, parameter :: pfb_waven_x = 2.
  real, parameter :: pfb_waven_y = 1.
  real, parameter :: pfb_waven_z = 0.
  real, parameter :: pfb_alpha_x = 0.


  !----------------------------------------------------------------------
  blkLimits = block%limits
  blkLimitsGC = block%limitsGC
  allocate(xLeft(blkLimitsGC(LOW, IAXIS):blkLimitsGC(HIGH, IAXIS)))
  allocate(xRight(blkLimitsGC(LOW, IAXIS):blkLimitsGC(HIGH, IAXIS)))
  allocate(xCenter(blkLimitsGC(LOW, IAXIS):blkLimitsGC(HIGH, IAXIS)))
  allocate(yCoord(blkLimitsGC(LOW, JAXIS):blkLimitsGC(HIGH, JAXIS)))
  allocate(zCoord(blkLimitsGC(LOW, KAXIS):blkLimitsGC(HIGH, KAXIS)))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  sizeX = SIZE(xLeft)
  sizeY = SIZE(yCoord)
  sizeZ = SIZE(zCoord)
  
  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, block, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, block, CENTER, gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, block, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, block, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, block, RIGHT_EDGE, gcell, xRight, sizeX)

#ifdef DEBUG_SIMULATION
98 format('initBlock:',A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format('initBlock:',A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print 99,"solnData" ,(lbound(solnData ,i),ubound(solnData ,i),i=1,4)
  print*,'blkLim  :',blkLimits
  print*,'blkLimGC:',blkLimitsGC
#endif
!------------------------------------------------------------------------------

  Lx = sim_xMax - sim_xMin
  Ly = sim_yMax - sim_yMin  
  Lz = sim_zMax - sim_zMin


  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(block%level,del)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(block,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Point to Blocks centered variables:
!  call Grid_getBlkPtr(blockID,solnData,CENTER)

!  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)


  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                real(i - NGUARD - 1)*del(IAXIS) + 0.5*del(IAXIS)

           yi = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                real(j - NGUARD - 1)*del(JAXIS) + 0.5*del(JAXIS)

           zi = coord(KAXIS) - 0.5*bsize(KAXIS) + &
                real(k - NGUARD - 1)*del(KAXIS) + 0.5*del(KAXIS)


           Phi_ijk = cos(2.*PI*xi*pfb_waven_x/Lx + pfb_alpha_x) * &
                     sin(2.*PI*yi*pfb_waven_y/Ly)*cos(2.*PI*zi*pfb_waven_z/Lz)

           
           F_ijk  = -4.*PI**2 * ( (pfb_waven_x/Lx)**2. + (pfb_waven_y/Ly)**2. + (pfb_waven_z/Lz)**2. ) * Phi_ijk
           
           solnData(ASOL_VAR,i,j,k) = Phi_ijk

           solnData(DENS_VAR,i,j,k) = F_ijk

        enddo
     enddo
  enddo


  ! set values for u,v velocities and pressure
  solnData(DIFF_VAR,:,:,:) = 0.0
  solnData(PFFT_VAR,:,:,:) = 0.0


!!$  write(*,*) 'BlockID=',blockID
!!$  write(*,*) 'Center coordinates=',coord
!!$  write(*,*) 'Size =',bsize
 

  ! Release pointer
!  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Simulation_initBlock
