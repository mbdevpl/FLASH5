!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( integer, intent(IN)    :: blockCount, 
!!               integer(:), intent(IN) :: blockList, 
!!               real, intent(IN)       ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
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

#include "Flash.h"

subroutine Burn (  dt  )    

  use bn_interface, ONLY : bn_burner   
  use bn_xnetData, ONLY : xnet_myid, xnet_nzbatchmx, xnet_inuc2unk
  use Burn_data, ONLY : bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_smallx, bn_useBurn, ionam
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getMaxRefinement, Grid_getLeafIterator, &
       Grid_releaseLeafIterator
  use Hydro_data, ONLY : hy_gcMaskSD
  use Hydro_interface, ONLY : Hydro_detectShock
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t

  !$ use omp_lib

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  !args
  real,    intent(in) :: dt

  ! locals
  integer :: blockID, thisBlock, blockCount
  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSize, jSize, kSize, iSizeGC, jSizeGC, kSizeGC
  integer :: iSize_max, jSize_max, kSize_max

  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
  logical, parameter :: getGuardCells = .true.

  real,    allocatable :: shock(:,:,:)
  real,    allocatable, target :: xIn(:,:,:,:,:), xOut(:,:,:,:,:)
  real,    allocatable, target :: sdot(:,:,:,:), tmp(:,:,:,:), rho(:,:,:,:)
  logical, allocatable, target :: burnedZone(:,:,:,:)

  real,    pointer :: xIn_batch(:,:,:), xOut_batch(:,:,:)
  real,    pointer :: sdot_batch(:,:), tmp_batch(:,:), rho_batch(:,:)
  logical, pointer :: burnedZone_batch(:,:)

  integer :: nzones, batchCount
  integer, dimension(:), allocatable :: sumBurn_TS_batch
  integer, dimension(:), allocatable :: batch_lo, batch_hi
  integer, dimension(:), allocatable :: sumBurn_TS

  real :: ei, ek, enuc
  integer :: i, j, k, m, n, ii, jj, kk, mm, nn

  integer :: level, maxLev
  type(leaf_iterator_t)  :: itor
  type(block_metadata_t) :: blockDesc

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

  !---------------------------------------------------------------------------------

  ! start the timer ticking
  call Timers_start("burn")

  call Timers_start("burn_top")

#ifdef FLASH_GRID_UG
     maxLev = 1
#else
     call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif

  blockCount = 0
  iSize_max = 1
  jSize_max = 1
  kSize_max = 1
  do level = 1, maxLev
     call Grid_getLeafIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)
        blockCount = blockCount+1
        blkLimits = blockDesc%limits
        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
        iSize_max = max(iSize_max,iSize)
        jSize_max = max(jSize_max,jSize)
        kSize_max = max(kSize_max,kSize)
        call itor%next()
     end do
     call Grid_releaseLeafIterator(itor)
  end do

  allocate(xIn(NSPECIES,iSize_max,jSize_max,kSize_max,blockCount))
  allocate(xOut(NSPECIES,iSize_max,jSize_max,kSize_max,blockCount))
  allocate(sdot(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(tmp(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(rho(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(burnedZone(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(batch_lo(blockCount))
  allocate(batch_hi(blockCount))
  allocate(sumBurn_TS(blockCount))

  burnedZone = .FALSE.

  if (.NOT. bn_useShockBurn) then
     call Grid_fillGuardCells(CENTER, ALLDIR, maskSize=NUNK_VARS, mask=hy_gcMaskSD)
  endif

  nzones = 0
  thisBlock = 0
  do level = 1, maxLev
     call Grid_getLeafIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)

        thisBlock = thisBlock + 1

        ! get dimensions/limits and coordinates
        blkLimitsGC = blockDesc%limitsGC
        iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        blkLimits = blockDesc%limits
        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

        !allocate(shock(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
        !               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
        !               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))

        ! identify the range of batches in each block (use floor/ceil in case of overlap)
        batch_lo(thisBlock) = nzones / xnet_nzbatchmx + 1
        nzones = nzones + iSize * jSize * kSize
        batch_hi(thisBlock) = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

        ! allocate space for dimensions
        allocate(xCoord(iSizeGC))
        allocate(yCoord(jSizeGC))
        allocate(zCoord(kSizeGC))

        call Grid_getCellCoords(IAXIS,blockDesc,CENTER,getGuardCells,xCoord,iSizeGC)
        call Grid_getCellCoords(JAXIS,blockDesc,CENTER,getGuardCells,yCoord,jSizeGC)
        call Grid_getCellCoords(KAXIS,blockDesc,CENTER,getGuardCells,zCoord,kSizeGC)

        ! Get a pointer to solution data 
        call Grid_getBlkPtr(blockDesc, solnData)

        ! Shock detector (post-poning AMREx implementation)
        !if (.NOT. bn_useShockBurn) then
        !   call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
        !      xCoord,yCoord,zCoord)
        !else
        !   shock(:,:,:) = 0
        !endif

        solnData(NMPI_VAR,:,:,:) = xnet_myid

        !$omp parallel do &
        !$omp collapse(3) &
        !$omp default(shared) &
        !$omp private(k,kk,j,jj,i,ii,okBurnTemp,okBurnDens,okBurnShock,okBurnNickel)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 kk = k - blkLimits(LOW,KAXIS) + 1
                 jj = j - blkLimits(LOW,JAXIS) + 1
                 ii = i - blkLimits(LOW,IAXIS) + 1

                 tmp(ii,jj,kk,thisBlock)  = solnData(TEMP_VAR,i,j,k)
                 rho(ii,jj,kk,thisBlock)  = solnData(DENS_VAR,i,j,k)
                 sdot(ii,jj,kk,thisBlock) = 0.0e0

                 ! Map the solution data into the order required by bn_burner
                 xIn(1:NSPECIES,ii,jj,kk,thisBlock) = solnData(xnet_inuc2unk,i,j,k)

                 okBurnTemp = .FALSE.
                 okBurnDens = .FALSE.
                 okBurnShock = .FALSE.
                 okBurnNickel = .FALSE.

                 okBurnTemp = (tmp(ii,jj,kk,thisBlock) >= bn_nuclearTempMin .AND. tmp(ii,jj,kk,thisBlock) <= bn_nuclearTempMax)
                 okBurnDens = (rho(ii,jj,kk,thisBlock) >= bn_nuclearDensMin .AND. rho(ii,jj,kk,thisBlock) <= bn_nuclearDensMax)
                 !okBurnShock = (shock(i,j,k) == 0.0 .OR. (shock(i,j,k) == 1.0 .AND. bn_useShockBurn))
                 okBurnShock = .TRUE.

                 if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

                    if (NI56_SPEC /= NONEXISTENT) then
                       okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                    else    ! nickel is not even a species in this simulation, so we'll always burn
                       okBurnNickel = .TRUE.
                    endif

                    if (okBurnNickel) then
                       burnedZone(ii,jj,kk,thisBlock) = .TRUE.
                    endif

                 endif

              end do
           end do
        end do
        !$omp end parallel do

        call Grid_releaseBlkPtr(blockDesc,solnData)
        nullify(solnData)

        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
        !deallocate(shock)

        call itor%next()

     end do
     call Grid_releaseLeafIterator(itor)
  end do

  call Timers_stop("burn_top")

  call Timers_start("burn_middle")

  ! get number of batches needed for all local zones (round up)
  batchCount = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

  ! reshape all local zone data arrays into batches
  tmp_batch (1:xnet_nzbatchmx,1:batchCount) => tmp (:,:,:,:)
  rho_batch (1:xnet_nzbatchmx,1:batchCount) => rho (:,:,:,:)
  sdot_batch(1:xnet_nzbatchmx,1:batchCount) => sdot(:,:,:,:)
  xIn_batch (1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xIn (:,:,:,:,:)
  xOut_batch(1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xOut(:,:,:,:,:)
  burnedZone_batch(1:xnet_nzbatchmx,1:batchCount) => burnedZone(:,:,:,:)

  allocate(sumBurn_TS_batch(batchCount))

  !$omp parallel do &
  !$omp schedule(runtime) &
  !$omp default(shared)
  do m = 1, batchCount
     ! Do the actual burn
     call bn_burner(dt, tmp_batch(:,m), rho_batch(:,m), xIn_batch(:,:,m), &
          xOut_batch(:,:,m), sdot_batch(:,m), burnedZone_batch(:,m), sumBurn_TS_batch(m))
  end do
  !$omp end parallel do

  ! get maximum sumBurn_TS_batch per block
  sumBurn_TS = 0
  do thisBlock = 1, blockCount
     sumBurn_TS(thisBlock) = maxval( sumBurn_TS_batch(batch_lo(thisBlock):batch_hi(thisBlock)) )
  end do

  deallocate(sumBurn_TS_batch)

  call Timers_stop("burn_middle")

  call Timers_start("burn_bottom")

  thisBlock = 0
  do level = 1, maxLev
     call Grid_getLeafIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)

        thisBlock = thisBlock + 1

        ! get dimensions/limits and coordinates
        blkLimits = blockDesc%limits

        ! Get a pointer to solution data 
        call Grid_getBlkPtr(blockDesc,solnData)

        ! Now put updated local data arrays back into unk through solnData pointer
        !$omp parallel do &
        !$omp collapse(3) &
        !$omp default(shared) &
        !$omp private(k,kk,j,jj,i,ii,ei,ek,enuc)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 kk = k - blkLimits(LOW,KAXIS) + 1
                 jj = j - blkLimits(LOW,JAXIS) + 1
                 ii = i - blkLimits(LOW,IAXIS) + 1

                 ! Map the solution data into the order required by bn_burner
                 solnData(xnet_inuc2unk,i,j,k) = xOut(1:NSPECIES,ii,jj,kk,thisBlock)

                 !  NOTE should probably do something here with eintSwitch for consistency
                 !  LBR will settle for simply using internal energy!
                 ! kinetic energy
                 ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                    solnData(VELY_VAR,i,j,k)**2 +  &
                    solnData(VELZ_VAR,i,j,k)**2)

                 ! internal energy, add on nuclear rate*timestep
                 enuc = dt*sdot(ii,jj,kk,thisBlock)
                 ei = solnData(ENER_VAR,i,j,k) + enuc - ek

#ifdef EINT_VAR
                 solnData(EINT_VAR,i,j,k) = ei
#endif
                 solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                 solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + enuc
#endif
                 solnData(ENUC_VAR,i,j,k) = sdot(ii,jj,kk,thisBlock)

              end do
           end do
        end do
        !$omp end parallel do
   
#ifdef FLASH_GRID_PARAMESH
        !bflags(1,blockID) = sumBurn_TS(thisBlock)
#endif
        solnData(MTSB_VAR,:,:,:) = sumBurn_TS(thisBlock)

        ! we've altered the EI, let's equilabrate
        call Timers_start("eos")
        if (any(burnedZone(:,:,:,thisBlock))) then

#ifdef FLASH_UHD_3T
           call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,solnData,CENTER) ! modified for 3T
#else 
           call Eos_wrapped(MODE_DENS_EI,blkLimits,solnData,CENTER)
#endif

        end if
        call Timers_stop("eos")

        call Grid_releaseBlkPtr(blockDesc,solnData)
        nullify(solnData)

        call itor%next()

     end do
     call Grid_releaseLeafIterator(itor)
  end do

  call Timers_stop("burn_bottom")

  deallocate(xIn)
  deallocate(xOut)
  deallocate(sdot)
  deallocate(tmp)
  deallocate(rho)
  deallocate(burnedZone)
  deallocate(batch_lo)
  deallocate(batch_hi)
  deallocate(sumBurn_TS)

  call Timers_stop("burn")

  return
end subroutine Burn
