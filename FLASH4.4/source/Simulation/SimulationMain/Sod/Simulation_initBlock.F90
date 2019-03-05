!!****if* source/Simulation/SimulationMain/Sod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!
!!***

!!REORDER(4): solnData
subroutine Simulation_initBlock(solnData, tileDesc)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos, &    
     &  sim_rhoLeft,  sim_pLeft, sim_uLeft, sim_rhoRight, sim_pRight, sim_uRight, &
     &  sim_smallX, sim_gamma, sim_smallP

#ifdef FLASH_3T
  use Simulation_data, ONLY : &
       sim_pionLeft, sim_peleLeft, sim_pradLeft, &
       sim_pionRight, sim_peleRight, sim_pradRight, &
       sim_gammaIon, sim_gammaEle
#endif
     
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax

  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone, gameZone, gamcZone

#ifdef SIMULATION_TWO_MATERIALS
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: mfrac
#endif

#ifdef FLASH_3T
  real :: peleZone, eeleZone
  real :: pionZone, eionZone
  real :: pradZone, eradZone
#endif

  integer :: lo(1:MDIM)
  integer :: hi(1:MDIM)

  lo(:) = tileDesc%limits(LOW,  :)
  hi(:) = tileDesc%limits(HIGH, :)
  allocate(  xLeft(lo(IAXIS):hi(IAXIS)))
  allocate( xRight(lo(IAXIS):hi(IAXIS)))
  allocate(xCenter(lo(IAXIS):hi(IAXIS)))
  allocate( yCoord(lo(JAXIS):hi(JAXIS)))
  allocate( zCoord(lo(KAXIS):hi(KAXIS)))
  xLeft = 0.0
  xRight = 0.0
  xCenter = 0.0
  yCoord = 0.0
  zCoord = 0.0

#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          lo, hi, zCoord)
#endif
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          lo, hi, yCoord)
#endif

  call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level, &
                          lo, hi, xLeft)
  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          lo, hi, xCenter)
  call Grid_getCellCoords(IAXIS, RIGHT_EDGE, tileDesc%level, &
                          lo, hi, xRight)

#ifdef DEBUG_SIMULATION
98 format('initBlock:',A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format('initBlock:',A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print 99,"solnData" ,(lbound(solnData ,i),ubound(solnData ,i),i=1,4)
  print*,'tile limits:',tileDesc%limits
  print*,'grown tile limits:',tileDesc%limitsGC
#endif
!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

  associate(lo => tileDesc%limits(LOW,  :), &
            hi => tileDesc%limits(HIGH, :))
  do k = lo(KAXIS), hi(KAXIS)
  
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
     
     do j = lo(JAXIS), hi(JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)

        ! The position of the shock in the current yz-row.
        lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        
        do i = lo(IAXIS), hi(IAXIS)
           
           ! get the cell center, left, and right positions in x
           xxL = xLeft(i)
           xx  = xCenter(i)
           xxR = xRight(i)
           
           ! initialize cells to the left of the initial shock.
           if (xxR <= lPosn) then
#ifdef FLASH_3T
              peleZone = sim_peleLeft
              pionZone = sim_pionLeft
              pradZone = sim_pradLeft
#else
              presZone = sim_pLeft
#endif              

              rhoZone = sim_rhoLeft
              velxZone = sim_uLeft * sim_xCos
              velyZone = sim_uLeft * sim_yCos
              velzZone = sim_uLeft * sim_zCos 

#ifdef SIMULATION_TWO_MATERIALS
              mfrac(LEFT_SPEC - SPECIES_BEGIN+1) = 1.0
              mfrac(RGHT_SPEC - SPECIES_BEGIN+1) = 0.0
#endif
              
              ! initialize cells which straddle the shock.  Treat them as though 1/2 of 
              ! the cell lay to the left and 1/2 lay to the right.
           elseif ((xxL < lPosn) .and. (xxR > lPosn)) then              
#ifdef FLASH_3T
              peleZone = 0.5 * (sim_peleLeft + sim_peleRight)
              pionZone = 0.5 * (sim_pionLeft + sim_pionRight)
              pradZone = 0.5 * (sim_pradLeft + sim_pradRight)
#else
              presZone = 0.5 * (sim_pLeft+sim_pRight)
#endif              
              
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoRight)
              velxZone = 0.5 *(sim_uLeft+sim_uRight) * sim_xCos
              velyZone = 0.5 *(sim_uLeft+sim_uRight) * sim_yCos
              velzZone = 0.5 *(sim_uLeft+sim_uRight) * sim_zCos
              
#ifdef SIMULATION_TWO_MATERIALS
              mfrac(LEFT_SPEC - SPECIES_BEGIN+1) = 0.5
              mfrac(RGHT_SPEC - SPECIES_BEGIN+1) = 0.5
#endif

              ! initialize cells to the right of the initial shock.
           else              
#ifdef FLASH_3T
              peleZone = sim_peleRight
              pionZone = sim_pionRight
              pradZone = sim_pradRight
#else
              presZone = sim_pRight
#endif              
              
              rhoZone = sim_rhoRight
              velxZone = sim_uRight * sim_xCos
              velyZone = sim_uRight * sim_yCos
              velzZone = sim_uRight * sim_zCos

#ifdef SIMULATION_TWO_MATERIALS
              mfrac(LEFT_SPEC - SPECIES_BEGIN+1) = 0.0
              mfrac(RGHT_SPEC - SPECIES_BEGIN+1) = 1.0
#endif
              
           endif

#ifdef FLASH_3T
           presZone = peleZone + pionZone + pradZone
#endif

           !put in default mass fraction values of all species
           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX

              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 solnData(n,i,j,k)= sim_smallX
              enddo
           end if

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
#ifdef SIMULATION_TWO_MATERIALS
           eosData(EOS_DENS) = rhoZone
           eosData(EOS_PRES) = presZone
           eosData(EOS_TEMP) = 1.0e8
           call Eos(MODE_DENS_PRES, 1, eosData, mfrac)
           eintZone = eosData(EOS_EINT)
           gameZone = 1.0+eosData(EOS_PRES)/eosData(EOS_DENS)/eosData(EOS_EINT)
           gamcZone = eosData(EOS_GAMC)
#else
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           gameZone = sim_gamma
           gamcZone = sim_gamma
#endif
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)

           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           

           solnData(DENS_VAR,i,j,k)=rhoZone
           solnData(DENS_VAR, i,j,k) =  rhoZone
           solnData(PRES_VAR, i,j,k) =  presZone
           solnData(VELX_VAR, i,j,k) =  velxZone
           solnData(VELY_VAR, i,j,k) =  velyZone
           solnData(VELZ_VAR, i,j,k) =  velzZone

#ifdef ENER_VAR
           solnData(ENER_VAR, i,j,k) =  enerZone
#endif
#ifdef EINT_VAR
           solnData(EINT_VAR, i,j,k) =  eintZone
#endif
#ifdef GAME_VAR          
           solnData(GAME_VAR, i,j,k) =  gameZone
#endif
#ifdef GAMC_VAR
           solnData(GAMC_VAR, i,j,k) =  gamcZone
#endif
#ifdef TEMP_VAR
# ifdef SIMULATION_TWO_MATERIALS
           solnData(TEMP_VAR, i,j,k) =  eosData(EOS_TEMP)
# else
           solnData(TEMP_VAR, i,j,k) =  1.e-10
# endif
#endif

#ifdef SIMULATION_TWO_MATERIALS
           solnData(LEFT_SPEC,i,j,k) = mfrac(LEFT_SPEC-SPECIES_BEGIN+1) 
           solnData(RGHT_SPEC,i,j,k) = mfrac(RGHT_SPEC-SPECIES_BEGIN+1) 
#endif

#ifdef FLASH_3T
           ! We must now compute the internal energy from the pressure
           ! for the ions, electrons, and radiation field:
           
           ! Electrons...
           eeleZone = peleZone / (sim_gammaEle - 1.0) / rhoZone
           eionZone = pionZone / (sim_gammaIon - 1.0) / rhoZone
           eradZone = 3.0 * pradZone / rhoZone
           
           solnData(EELE_VAR, i,j,k) =  eeleZone
           solnData(EION_VAR, i,j,k) =  eionZone
           solnData(ERAD_VAR, i,j,k) =  eradZone
#ifdef DFCF_VAR
           solnData(DFCF_VAR, i,j,k) =  eintZone
#endif
           eintZone = eeleZone + eionZone + eradZone !recompute
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP/rhoZone)
#ifdef ENER_VAR
           solnData(ENER_VAR, i,j,k) =  enerZone
#endif
#ifdef EINT_VAR
           solnData(EINT_VAR, i,j,k) =  eintZone
#endif
#ifdef GAME_VAR
           solnData(GAME_VAR, i,j,k) =  presZone/(rhoZone*eintZone)+1.0
#endif
#ifdef PION_VAR
           solnData(PION_VAR, i,j,k) =  pionZone
#endif
#ifdef PELE_VAR
           solnData(PELE_VAR, i,j,k) =  peleZone
#endif
#ifdef PRAD_VAR
           solnData(PRAD_VAR, i,j,k) =  pradZone
#endif
#endif
        enddo
     enddo
  enddo
  end associate

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock
