!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( real, intent(IN) ::  dt  )
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   dt  --       passed to the internal bn_burner module
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useBurnTable -- Boolean, False.  Controls the generation of reaction rates.
!!                TRUE interpolates from a stored table; FALSE generates them
!!                analytically.
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  algebra -- Integer, 1, [1,2].  Controls choice of linear algebra package used
!!                for matrix solution.  1=Ma28 sparse package, 2=Gift hardwired package.
!!  odeStepper -- Integer, 1, [1,2].  Controls time integration routines.
!!                1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!  nuclearTempMin/Max -- Real, 1.1E+8/1.0E+12.  Minimum and maximum temperature
!!                ranges where burning can occur
!!  nuclearDensMin/Max -- Real, 1.0E-10/1.0E+14.  Minimum and maximum density range
!!                where burning can occur.
!!  nuclearNI56Max -- Real, 1.0.  Maximum mass fraction of nickel where burning
!!                can occur.
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.See Burn_computeDt for details.
!!
!! NOTES
!!
!!  The burning unit adds a new mesh variable ENUC_VAR which is the nuclear energy
!!             generation rate
!!
!!***

!!REORDER(4): solnData

subroutine Burn (  dt  )

  use Burn_data, ONLY : bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_smallx, bn_useBurn
  use bn_interface, ONLY : bn_mapNetworkToSpecies, bn_burner

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getCellCoords, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getMaxRefinement, &
       Grid_getLeafIterator, Grid_releaseLeafIterator
  use Eos_interface, ONLY : Eos_wrapped
  use Hydro_data, ONLY : hy_gcMaskSD
  use Hydro_interface, ONLY : Hydro_shockStrength

  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t

  !$ use omp_lib

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  !args
  real, intent(in) :: dt

  ! locals
  integer :: i, j, k, n, speciesMap

  real, dimension(NSPECIES) :: xIn, xOut
  real :: sdot, tmp, rho, ei, ek, enuc

  logical :: burnedZone
  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
  logical, parameter :: getGuardCells = .true.

  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSize, jSize, kSize, iSizeGC, jSizeGC, kSizeGC

  integer, parameter :: shock_mode = 1
  real, parameter :: shock_thresh = 0.33
  real, allocatable, dimension(:,:,:) :: shock

  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: level, maxLev
  type(leaf_iterator_t)  :: itor
  type(block_metadata_t) :: blockDesc

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

  !---------------------------------------------------------------------------------

  ! start the timer ticking
  call Timers_start("burn")

#ifdef FLASH_GRID_UG
     maxLev = 1
#else
     call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif

  if (.NOT. bn_useShockBurn) then
     call Grid_fillGuardCells(CENTER, ALLDIR, maskSize=NUNK_VARS, mask=hy_gcMaskSD)
  endif

  do level = 1, maxLev
     call Grid_getLeafIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)

        burnedZone = .FALSE.

        ! get dimensions/limits and coordinates
        blkLimitsGC = blockDesc%limitsGC
        iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        blkLimits = blockDesc%limits
        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

        ! allocate space for dimensions
        allocate(xCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
        allocate(yCoord(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
        allocate(zCoord(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

        allocate(shock(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

        call Grid_getCellCoords(IAXIS,blockDesc,CENTER,getGuardCells,xCoord,iSizeGC)
        call Grid_getCellCoords(JAXIS,blockDesc,CENTER,getGuardCells,yCoord,jSizeGC)
        call Grid_getCellCoords(KAXIS,blockDesc,CENTER,getGuardCells,zCoord,kSizeGC)

        ! Get a pointer to solution data
        call Grid_getBlkPtr(blockDesc, solnData)

        ! Shock detector
        if (.NOT. bn_useShockBurn) then
           call Hydro_shockStrength(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
              xCoord,yCoord,zCoord,shock_thresh,shock_mode)
        else
           shock(:,:,:) = 0.0
        endif

        !$omp parallel do &
        !$omp collapse(3) &
        !$omp default(shared) &
        !$omp private(i,j,k,n,speciesMap,tmp,rho,sdot,xIn,xOut,ei,ek,enuc, &
        !$omp         okBurnTemp,okBurnDens,okBurnShock,okBurnNickel)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 tmp  = solnData(TEMP_VAR,i,j,k)
                 rho  = solnData(DENS_VAR,i,j,k)
                 sdot = 0.0e0

                 okBurnTemp = (tmp >= bn_nuclearTempMin .AND. tmp <= bn_nuclearTempMax)
                 okBurnDens = (rho >= bn_nuclearDensMin .AND. rho <= bn_nuclearDensMax)
                 okBurnShock = (shock(i,j,k) <= 0.0 .OR. (shock(i,j,k) > 0.0 .AND. bn_useShockBurn))

                 if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

                    if (NI56_SPEC /= NONEXISTENT) then
                       okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                    else    ! nickel is not even a species in this simulation, so we'll always burn
                       okBurnNickel = .TRUE.
                    endif

                    if (okBurnNickel) then

                       !$omp atomic write
                       burnedZone = .TRUE.
                       !$omp end atomic

                       ! Map the solution data into the order required by bn_burner
                       do n = 1, NSPECIES
                          call bn_mapNetworkToSpecies(n,speciesMap)
                          xIn(n) = solnData(speciesMap,i,j,k)
                       end do

                       ! Do the actual burn
                       call bn_burner(dt, tmp, rho, xIn, xOut, sdot)

                       !  Map it back NOTE someday make a nicer interface....
                       do n = 1, NSPECIES
                          call bn_mapNetworkToSpecies(n,speciesMap)
                          solnData(speciesMap,i,j,k) = xOut(n)
                       end do

                       ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                            solnData(VELY_VAR,i,j,k)**2 +  &
                            solnData(VELZ_VAR,i,j,k)**2)

                       ! internal energy, add on nuclear rate*timestep
                       enuc = dt*sdot
                       ei = solnData(ENER_VAR,i,j,k) + enuc - ek

#ifdef EINT_VAR
                       solnData(EINT_VAR,i,j,k) = ei
#endif
                       solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                       solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + dt*sdot
#endif
                       solnData(ENUC_VAR,i,j,k) = sdot

                    endif
                 endif

              end do
           end do
        end do
        !$omp end parallel do

        ! we've altered the EI, let's equilabrate
        if (burnedZone) then

#ifdef FLASH_UHD_3T
           call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,solnData,CENTER) ! modified for 3T
#else
           call Eos_wrapped(MODE_DENS_EI,blkLimits,solnData,CENTER)
#endif

        end if

        call Grid_releaseBlkPtr(blockDesc,solnData)
        nullify(solnData)

        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
        deallocate(shock)

        call itor%next()

     end do
     call Grid_releaseLeafIterator(itor)
  end do

  call Timers_stop("burn")

  return
end subroutine Burn
