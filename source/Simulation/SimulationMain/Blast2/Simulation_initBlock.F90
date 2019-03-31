!!****if* source/Simulation/SimulationMain/Blast2/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!                       
!!
!!
!! DESCRIPTION
!!  Initial conditions for setting up the Woodward Collela 
!!  two blast-wave problem.
!!  Reference:   Woodward, P. & Colella, P. 1984, J. Comp. Phys., 54, 115
!!
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!!
!!***


subroutine Simulation_initBlock(solnData,blockDesc)

  use Simulation_data, ONLY: sim_xCos, sim_yCos, sim_zCos, &
     & sim_posnL, sim_posnR, sim_smallX, sim_gamma, sim_smallP, &
       sim_smallE, sim_smallRho, &
     & sim_rhoLeft, sim_rhoMid, sim_rhoRight, &
     & sim_pLeft, sim_pMid, sim_pRight, &
     & sim_uLeft, sim_uMid, sim_uRight, &
     sim_eMassInUAmu
  use Grid_interface, ONLY :Grid_getCellCoords
  use block_metadata, ONLY : block_metadata_t

#include "constants.h"
#include "Flash.h"


  implicit none      


  real,                   pointer    :: solnData(:,:,:,:)
  type(block_metadata_t), intent(in) :: blockDesc

! compute the maximum length of a vector in each coordinate direction 
! (including guardcells)

  integer :: i, j, k, n
  
  real :: xx, yy, zz, xxL, xxR
  real :: lPosn0L, lPosnL, lPosn0R, lPosnR
  real :: rhoZone, velXZone, velYZone, velzZone, presZone, & 
          eintZone, enerZone, ekinZone
  real :: Zfree, relA, YeFree
  real :: ionEnerFrac, eleEnerFrac
  real :: ionMassFrac, eleMassFrac, radMassFrac
  
  real, allocatable,dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  logical :: gcell = .true.

  blkLimitsGC(:,:) = blockDesc%limitsGC
  blkLimits(:,:) = blockDesc%limits
  
  allocate(xLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xRight(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xCenter(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(yCoord(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(zCoord(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0
  sizeX = SIZE(xCenter)
  sizeY = SIZE(yCoord)
  sizeZ = SIZE(zCoord)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockDesc, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockDesc, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockDesc, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockDesc, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockDesc, RIGHT_EDGE, gcell, xRight, sizeX)


! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

     zz = zCoord(k) ! get cell center coordinates in  z-direction
     axis(KAXIS) = k
     ! Where along the x-axis the shock intersects the xz-plane at the current z
     lPosn0L = sim_posnL - zz*sim_zCos/sim_xCos
     lPosn0R = sim_posnR - zz*sim_zCos/sim_xCos
     
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS) = j

        yy = yCoord(j)! get cell center coordinates in the y-direction

        ! The position of the shock in the current yz-row.

        lPosnL = lPosn0L - yy*sim_yCos/sim_xCos
        lPosnR = lPosn0R - yy*sim_yCos/sim_xCos


        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS) = i
           
           xx  = xCenter(i) ! get cell center, left, and right positions in x
           
           xxL = xLeft(i)
           xxR = xRight(i)
           
           ! Initialize cells in the left part of the grid.
           
           if (xxR < lPosnL) then
              
              rhoZone = sim_rhoLeft
              presZone = sim_pLeft
              
              velXZone = sim_uLeft * sim_xCos
              velYZone = sim_uLeft * sim_yCos
              velzZone = sim_uLeft * sim_zCos
              
              ! Initialize cells which straddle the left shock. 
              ! Treat them as though 
              ! 1/2 of the cell lay to the left and 1/2 lay to the right.
              
           elseif ((xxL < lPosnL) .and. (xxR > lPosnL)) then
              
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoMid)
              presZone = 0.5 * (sim_pLeft+sim_pMid)
              
              velXZone = 0.5 * (sim_uLeft+sim_uMid) * sim_xCos
              velYZone = 0.5 * (sim_uLeft+sim_uMid) * sim_yCos
              velzZone = 0.5 * (sim_uLeft+sim_uMid) * sim_zCos
              
              ! Initialize cells in the middle part of the grid.
              
           elseif ((xxL > lPosnL) .and. (xxR < lPosnR)) then
              
              rhoZone = sim_rhoMid
              presZone = sim_pMid
              
              velXZone = sim_uMid * sim_xCos
              velYZone = sim_uMid * sim_yCos
              velzZone = sim_uMid * sim_zCos
              
              ! Initialize cells which straddle the right shock.
              
           elseif ((xxL < lPosnR) .and. (xxR > lPosnR)) then
              
              rhoZone = 0.5 * (sim_rhoRight+sim_rhoMid)
              presZone = 0.5 * (sim_pRight+sim_pMid)
              
              velXZone = 0.5 * (sim_uRight+sim_uMid) * sim_xCos
              velYZone = 0.5 * (sim_uRight+sim_uMid) * sim_yCos
              velzZone = 0.5 * (sim_uRight+sim_uMid) * sim_zCos
              
              ! Initialize cells in the right part of the grid.
              
           else
              
              rhoZone = sim_rhoRight
              presZone = sim_pRight
              
              velXZone = sim_uRight * sim_xCos
              velYZone = sim_uRight * sim_yCos
              velzZone = sim_uRight * sim_zCos
              
           endif
           

#if NSPECIES > 0
           !put default species
           solnData(SPECIES_END,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX

           !if there is only 1 species this loop will not execute
           do n = SPECIES_BEGIN, SPECIES_END-1
              solnData(n,i,j,k)=sim_smallX
           enddo
#endif
           
           ekinZone = 0.5 * (velXZone**2 + velYZone**2 + velzZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           Zfree = 6
           relA = 12
           YeFree = 6/relA
           ionEnerFrac = 1.0   / (Zfree + 1.0)
           eleEnerFrac = Zfree / (Zfree + 1.0)

           eleMassFrac = Zfree*sim_eMassInUAmu / relA
           ionMassFrac = 1.0 - eleMassFrac

           radMassFrac = 0.0

           solnData(DENS_VAR,i,j,k)=rhoZone
           solnData(PRES_VAR,i,j,k)=presZone
           solnData(VELX_VAR,i,j,k)=velxZone
           solnData(VELY_VAR,i,j,k)=velyZone
           solnData(VELZ_VAR,i,j,k)=velzZone
#ifdef GAME_VAR          
           solnData(GAME_VAR,i,j,k)=sim_gamma
#endif

#ifdef GAMC_VAR
           solnData(GAMC_VAR,i,j,k)=sim_gamma
#endif

#ifdef ENER_VAR
           solnData(ENER_VAR,i,j,k)=enerZone
#endif
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k)=eintZone
#endif
#ifdef TEMP_VAR
           solnData(TEMP_VAR,i,j,k)= 1e-6
#endif

           
#ifdef PION_VAR
           solnData(PION_VAR,i,j,k)=presZone*ionEnerFrac
#endif
#ifdef PELE_VAR
           solnData(PELE_VAR,i,j,k)=presZone*elecEnerFrac
#endif
#ifdef E1_VAR
           solnData(E1_VAR,i,j,k)= ekinZone*ionMassFrac+eintZone*ionEnerFrac
#endif
#ifdef E2_VAR
           solnData(E2_VAR,i,j,k)= ekinZone*eleMassFrac+eintZone*eleEnerFrac
#endif
#ifdef E3_VAR
           solnData(E3_VAR,i,j,k)=sim_smallE
#endif
#ifdef EION_VAR
           solnData(EION_VAR,i,j,k)=eintZone*ionEnerFrac
#endif
#ifdef EELE_VAR
           solnData(EELE_VAR,i,j,k)=eintZone*eleEnerFrac
#endif
#ifdef ERAD_VAR
           solnData(ERAD_VAR,i,j,k)=sim_smallE
#endif
#ifdef TION_VAR
           solnData(TION_VAR,i,j,k)= 1e-6
#endif
#ifdef TELE_VAR
           solnData(TELE_VAR,i,j,k)= 1e-6
#endif
#ifdef TRAD_VAR
           solnData(TRAD_VAR,i,j,k)= 1e-6
#endif
#ifdef YE_MSCALAR
!!$           solnData(YE_MSCALAR,i,j,k)=0.5
           solnData(YE_MSCALAR,i,j,k)=YeFree
#endif
#ifdef SUMY_MSCALAR
!!$           solnData(SUMY_MSCALAR,i,j,k)= 7.292E-02
           solnData(SUMY_MSCALAR,i,j,k)=1.0/relA
#endif

        enddo
     enddo
  enddo
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)
  
  return
end subroutine Simulation_initBlock

