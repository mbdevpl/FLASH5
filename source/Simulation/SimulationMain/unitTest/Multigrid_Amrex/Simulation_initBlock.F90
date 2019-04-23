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
!!  a specified tile.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  tile -          the tile to update
!!  
!!
!! 
!!
!!***

!!-----!! Do not REORDER(4): solnData

subroutine Simulation_initBlock(solnData,tileDesc)

!  use Simulation_data
  use Simulation_data, ONLY :sim_xMin,sim_xMax,sim_yMin,sim_yMax,sim_zMin,sim_zMax
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc
  integer :: tileDescID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(MDIM) :: lo, hi

  real,allocatable, dimension(:) ::xCenter,yCenter,zCenter
  integer :: sizeX,sizeY,sizeZ

  real :: Lx, Ly, Lz, xi, yi, zi, Phi_ijk, F_ijk


  real :: pfb_waven_x
  real :: pfb_waven_y
  real :: pfb_waven_z
  real :: pfb_alpha_x

  logical :: gcell = .true.

  !----------------------------------------------------------------------
  lo=tileDesc%limits(LOW,:)
  hi=tileDesc%limits(HIGH,:)
  pfb_waven_x = 2.
  pfb_waven_y = 1.
  pfb_waven_z = 2.
  pfb_alpha_x = 0.
  if (NDIM < 3) pfb_waven_z=0.
  if (NDIM < 2) pfb_waven_y=0.
  allocate(xCenter(lo(IAXIS):hi(IAXIS)))
  allocate(yCenter(lo(JAXIS):hi(JAXIS)))
  allocate(zCenter(lo(KAXIS):hi(KAXIS)))
  xCenter = 0.0
  yCenter = 0.0
  zCenter = 0.0

  sizeX = SIZE(xCenter)
  sizeY = SIZE(yCenter)
  sizeZ = SIZE(zCenter)
  
  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)
  if (NDIM == 3) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

#ifdef DEBUG_SIMULATION
98 format('initTile:',A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format('initTile:',A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print 99,"solnData" ,(lbound(solnData ,i),ubound(solnData ,i),i=1,4)
  print*,'blkLim  :',lo,hi
#endif
!------------------------------------------------------------------------------

  Lx = sim_xMax - sim_xMin
  Ly = sim_yMax - sim_yMin  
  Lz = sim_zMax - sim_zMin

  do       k = lo(KAXIS), hi(KAXIS)
     do    j = lo(JAXIS), hi(JAXIS)
        do i = lo(IAXIS), hi(IAXIS)
          xi=xCenter(i)
          yi=yCenter(j)
          zi=zCenter(k)

           Phi_ijk = cos(2.*PI*xi*pfb_waven_x/Lx + pfb_alpha_x) * &
                     sin(2.*PI*yi*pfb_waven_y/Ly)*cos(2.*PI*zi*pfb_waven_z/Lz)

  
           F_ijk  = -4.*PI**2 * ( (pfb_waven_x/Lx)**2. + (pfb_waven_y/Ly)**2. + (pfb_waven_z/Lz)**2. ) * Phi_ijk
           
           solnData(i,j,k,ASOL_VAR) = Phi_ijk

           solnData(i,j,k,RHS_VAR) = F_ijk

        enddo
     enddo
  enddo


  ! set values for other variables
  solnData(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS), DIFF_VAR) = 0.0
  solnData(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS), NSOL_VAR) = 0.0

  return

end subroutine Simulation_initBlock
