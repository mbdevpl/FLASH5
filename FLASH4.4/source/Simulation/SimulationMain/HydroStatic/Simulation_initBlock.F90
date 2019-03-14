!!****if* source/Simulation/SimulationMain/HydroStatic/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer     :: solnData(:,:,:,:),
!!                            Grid_tile_t(IN)  :: tileDesc  )
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Hydrostatic boundary 
!!  condition test problem.
!!
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  tileDesc -        describes the tile or block to initialize
!!
!!
!!***

!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_SIMULATION
#endif

subroutine Simulation_initBlock(solnData,tileDesc)
  ! get the needed unit scope data
  use Grid_interface, ONLY : Grid_getCellCoords,&
                              Grid_putPointData
  use Simulation_data, ONLY : sim_gamma, sim_xyzRef, sim_presRef, &
       sim_gravVector, sim_gravDirec, &
       sim_xyzRef, sim_densRef, sim_tempRef, &
       sim_molarMass, sim_gasconstant
  use Simulation_data, ONLY : myPE => sim_meshMe
  use Grid_tile, ONLY : Grid_tile_t 
  implicit none


! get all the constants
#include "constants.h"
#include "Flash.h"

! declare arguments and indicate whether they are input or output
  real,                   pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t),      intent(in) :: tileDesc




! declare all local variables.
  integer :: i, j, k, n
  real :: xx, yy,  zz, gh

  ! arrays to hold coordinate information for the block
  real,allocatable, dimension(:) ::xCoord,yCoord,zCoord


  ! array to get integer indices defining the beginning and the end
  ! of a tile (or block).
  integer, dimension(2,MDIM) :: tileLimits


  integer, dimension(MDIM) :: axis
  integer,parameter :: dataSize = 1  ! for Grid put data function
  logical,parameter :: gcell = .true.


  ! these variables store the calculated initial values of physical
  ! variables a grid point at a time.
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, &
       enerZone, ekinZone, tempZone

  ! get the integer endpoints of the block in all dimensions
  tileLimits = tileDesc%limits

  allocate(xCoord(tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS))); xCoord = 0.0
  allocate(yCoord(tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS))); yCoord = 0.0
  allocate(zCoord(tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS))); zCoord = 0.0

  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          tileLimits(LOW,  :), &
                          tileLimits(HIGH, :), &
                          xCoord)
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          tileLimits(LOW,  :), &
                          tileLimits(HIGH, :), &
                          yCoord)
#endif
#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          tileLimits(LOW,  :), &
                          tileLimits(HIGH, :), &
                          zCoord)
#endif


  !-----------------------------------------------------------------------------
  ! loop over all of the zones in the current block and set the variables.
  !-----------------------------------------------------------------------------
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     zz = zCoord(k) ! coordinates of the cell center in the z-direction

     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        yy = yCoord(j) ! center coordinates in the y-direction

        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           xx = xCoord(i)

           velxZone = 0.0
           velyZone = 0.0
           velzZone = 0.0

           axis(IAXIS) = i   ! Get the position of the cell in the block
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! Compute the gas energy and set the gamma-values
           ! needed for the equation of  state.

           tempZone = sim_tempRef

           gh = sim_gravVector(1) * (xx-sim_xyzRef) + &
                 sim_gravVector(2) * (yy-sim_xyzRef) + &
                 sim_gravVector(3) * (zz-sim_xyzRef)

!           print*,i,j,exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
           rhoZone = sim_densRef * exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
           presZone = sim_presRef * rhoZone / sim_densRef

           ekinZone = 0.0

           enerZone = presZone / (sim_gamma-1.)
           enerZone = enerZone / rhoZone
           enerZone = enerZone + ekinZone

           ! store the variables in the current zone via the Grid_putPointData method

           solnData(DENS_VAR,i,j,k) = rhoZone
           solnData(PRES_VAR,i,j,k) = presZone
           solnData(VELX_VAR,i,j,k) = velxZone
           solnData(VELY_VAR,i,j,k) = velyZone
           solnData(VELZ_VAR,i,j,k) = velzZone
           solnData(ENER_VAR,i,j,k) = enerZone
           solnData(EINT_VAR,i,j,k) = enerZone
           solnData(TEMP_VAR,i,j,k) = tempZone
           solnData(GAME_VAR,i,j,k) = sim_gamma
           solnData(GAMC_VAR,i,j,k) = sim_gamma

        end do
     end do
  end do

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)


end subroutine Simulation_initBlock
