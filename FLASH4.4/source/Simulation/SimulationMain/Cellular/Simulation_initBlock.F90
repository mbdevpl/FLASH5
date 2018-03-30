!!****if* source/Simulation/SimulationMain/Cellular/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer, INTENT(in)::blockID)
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a planar detonation,
!!   perturbed with some noise.  This detonation should be unstable
!!   to becoming cellular at the detonation front.
!!
!! ARGUMENTS
!!
!!  blockID:        integer  the number of the block to initialize
!!
!!
!! PARAMETERS
!!
!!    xhe4               mass fraction of he4
!!    xc12               mass fraction of c12
!!    xo16               mass fraction of o16
!!    rhoAmbient         density of the cold upstream material 
!!    tempAmbient        temperature of the cold upstream material
!!    velxAmbient        x-velocity of the cold upstream material
!!    rhoPerturb         density of the post shock material
!!    tempPerturb        temperature of the post shock material
!!    velxPerturb        x-velocity of the post shock material
!!    radiusPerturb      distance below which the perturbation is applied
!!    xCenterPerturb     origin of the of the perturbation
!!    yCenterPerturb     origin of the of the perturbation
!!    zCenterPerturb     origin of the of the perturbation
!!    usePseudo1d        .true. for a 1d initial configuration, with the ??
!!                          copied along the y and z directions
!!                       .false. for a spherical configuration
!!    noiseAmplitude     amplitude of the white noise added to the perturbation
!!    noiseDistance      distances above and below radiusPerturb get noise added
!!    xmax               boundary of domain
!!    xmin               boundary of domain
!!    ymax               boundary of domain
!!    ymin               boundary of domain
!!    zmax               boundary of domain
!!    zmin               boundary of domain
!!
!!    smallx             smallest allowed abundance
!!    smlrho             smallest allowed density
!!
!!  NOTES
!!    Species used in this simulation are HE4 (helium-4), C12 (carbon-12), O16 (oxygen-16)
!!
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!  
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(solnData,block)

  use Simulation_data, ONLY: sim_smallRho, sim_smallx, sim_radiusPerturb, sim_usePseudo1d, &
     sim_xhe4, sim_xc12, sim_xo16, &
     sim_rhoAmbient, sim_tempAmbient, sim_velxAmbient, &
     sim_rhoPerturb, sim_tempPerturb, sim_velxPerturb, &
     sim_noiseAmplitude, sim_noiseDistance, &
     sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb, &
     sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax  
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Grid_interface, ONLY : Grid_getCellCoords
  use Eos_interface, ONLY : Eos
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block

!  Local variables

  real :: xx, yy, zz, dist
  logical, parameter :: useGuardCell = .TRUE.

  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction

  real,allocatable,dimension(:) :: xCoordsCell,yCoordsCell,zCoordsCell
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ

  integer :: i, j, k, n

  integer :: icount
  integer, parameter :: ifail = -1
  real, allocatable  :: rvec(:)                   ! for the random number generator
  integer            :: rvecSize=0             ! number of random numbers needed,
                                                     ! calculated below
  integer, parameter :: iseed = -867690
  integer            :: iseedUse

  ! variables needed for the eos call
  real :: temp_zone, rho_zone, vel_zone
  real :: ptot, eint, etot
  real, dimension(EOS_NUM)  :: eosData

  iseedUse = iseed
! ----------------------------------------------------------------------------------------------

  ! Get the indices of the blocks
  blkLimits = block%limits
  blkLimitsGC = block%limitsGC
  allocate(xCoordsCell(blkLimitsGC(LOW, IAXIS):blkLimitsGC(HIGH, IAXIS))); xCoordsCell = 0.0
  allocate(yCoordsCell(blkLimitsGC(LOW, JAXIS):blkLimitsGC(HIGH, JAXIS))); yCoordsCell = 0.0
  allocate(zCoordsCell(blkLimitsGC(LOW, KAXIS):blkLimitsGC(HIGH, KAXIS))); zCoordsCell = 0.0
  sizeX = SIZE(xCoordsCell)
  sizeY = SIZE(yCoordsCell)
  sizeZ = SIZE(zCoordsCell)

  call Grid_getCellCoords(KAXIS, block, CENTER, useGuardCell, zCoordsCell, sizeZ)
  call Grid_getCellCoords(JAXIS, block, CENTER, useGuardCell, yCoordsCell, sizeY)
  call Grid_getCellCoords(IAXIS, block, CENTER, useGuardCell, xCoordsCell, sizeX)

  ! the initial composition
  massFraction(:)    = sim_smallx 
  if (HE4_SPEC > 0) massFraction(HE4_SPEC) = max(sim_xhe4,sim_smallx)
  if (C12_SPEC > 0) massFraction(C12_SPEC) = max(sim_xc12,sim_smallx)
  if (O16_SPEC > 0) massFraction(O16_SPEC) = max(sim_xo16,sim_smallx)

  !..get a blocks worth of random numbers between 0.0 and 1.0
  rvecSize = sizeX*sizeY*sizeZ
  allocate(rvec(rvecSize))
  call sim_ranmar(iseedUse, rvec, rvecSize)

  icount = 0

  ! now fill the master arrays

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     if (NDIM == 3) zz = zCoordsCell(k)

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        if (NDIM >= 2) yy = yCoordsCell(j)

        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx = xCoordsCell(i)
           icount = icount + 1

           ! compute the distance from the center
           if (NDIM == 1) then
              dist = xx - sim_xCenterPerturb
           else if (NDIM == 2) then
              if (sim_usePseudo1d) then
                 dist = xx - sim_xCenterPerturb
              else
                 dist = sqrt((xx - sim_xCenterPerturb)**2 + & 
                      &                 (yy - sim_yCenterPerturb)**2)
              endif
           else if (NDIM == 3) then
              if (sim_usePseudo1d) then
                 dist = xx - sim_xCenterPerturb
              else
                 dist = sqrt((xx - sim_xCenterPerturb)**2 + & 
                      &                 (yy - sim_yCenterPerturb)**2 + & 
                      &                 (zz - sim_zCenterPerturb)**2)
              endif
           endif

           ! set the temperature, density, and x-velocity
           if (dist <= sim_radiusPerturb) then
              temp_zone = sim_tempPerturb
              rho_zone  = sim_rhoPerturb
              vel_zone  = sim_velxPerturb
           else
              temp_zone = sim_tempAmbient
              rho_zone  = sim_rhoAmbient
              vel_zone  = sim_velxAmbient
           endif


           !..seed the initial conditions with some white noise
           if ( abs(dist - sim_radiusPerturb) <= sim_noiseDistance) then
            !!  print *, 'sim_noiseAmplitude = ',sim_noiseAmplitude
            !!  print *,'rvec, icount', icount, rvec(icount)
              rho_zone = rho_zone *  & 
                   &               (1.0 + sim_noiseAmplitude * (1.0 - 2.0 * rvec(icount)))
           end if


           !  Need input of density and temperature
           eosData(EOS_TEMP) = temp_zone
           eosData(EOS_DENS) = rho_zone

           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)

           temp_zone = eosData(EOS_TEMP)
           rho_zone = eosData(EOS_DENS)
           ptot = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)

           ! calculate kinetic energy and total energy
           !! this was NOT done in flash2
           etot = eint + 0.5*vel_zone**2

           ! store the values
           ! fill the flash arrays

           solnData(TEMP_VAR,i,j,k)=temp_zone
           solnData(DENS_VAR,i,j,k)=rho_zone
           solnData(PRES_VAR,i,j,k)=ptot
           solnData(EINT_VAR,i,j,k)=e
           solnData(ENER_VAR,i,j,k)=etot
           solnData(GAMC_VAR,i,j,k)=eosData(EOS_GAMC)
           solnData(GAME_VAR,i,j,k)=(ptot/(etot*sim_rhoAmbient) + 1.0)
           solnData(VELX_VAR,i,j,k)=vel_zone
           solnData(VELY_VAR,i,j,k)=0.0
           solnData(VELZ_VAR,i,j,k)=0.0
           do n = SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k)=massFraction(n)
           enddo

           !..end of 3d loops
        enddo  ! end of k loop
     enddo     ! end of j loop
  enddo        ! end of i loop

  ! cleanup
  deallocate(rvec)
  deallocate(xCoordsCell)
  deallocate(yCoordsCell)
  deallocate(zCoordsCell)
  return
end subroutine Simulation_initBlock
