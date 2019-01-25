#include "Flash.h"
#include "constants.h"

subroutine hy_advance(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getDeltas, &
                                  Grid_zeroFluxData, &
                                  Grid_getTileIterator, &
                                  Grid_releaseTileIterator, &
                                  Grid_getMaxRefinement, Grid_conserveFluxes, Grid_putFluxData
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use hy_interface,        ONLY : hy_computeFluxes, hy_updateSolution
  use hy_memInterface,     ONLY : hy_memAllocScratch,      &
                                  hy_memDeallocScratch
  use flash_iterator,      ONLY : flash_iterator_t
  use flash_tile,          ONLY : flash_tile_t
  use Hydro_data,          ONLY : hy_fluxCorrect, hy_fluxCorrectPerLevel

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  integer,save :: sweepDummy = SWEEP_ALL

  integer:: level, maxLev
  
  type(flash_iterator_t)  :: itor
  type(flash_tile_t)      :: tileDesc

  nullify(Uin)
  nullify(Uout)

  call Grid_zeroFluxData

#ifdef USE_NOFLUXCORR_SHORTCUT
  ! ***** FIRST VARIANT: OPTIMIZED (somewhat) FOR hy_fluxCorrect==.FALSE. *****
  if (.NOT. hy_fluxCorrect) then
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     call Timers_start("hy_advance")
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        level            = tileDesc%level
        blkLimits(:,:)   = tileDesc%limits
        blkLimitsGC(:,:) = tileDesc%limitsGC

        call tileDesc%getDataPtr(Uin, CENTER)
        call tileDesc%deltas(del)
        Uout => Uin             ! hy_computeFluxes will ALSO update the solution through the Uout pointer!
        call hy_computeFluxes(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
!!$        call hy_updateSolution(blockDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uin, CENTER)
        nullify(Uout)
        call itor%next()
     end do
     call Timers_stop("hy_advance")
     call Grid_releaseTileIterator(itor)

     RETURN                     ! DONE, return from here!
  end if
#endif

  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.

!!$  call hy_memAllocScratch(SCRATCH_CTR,HY_VAR1_SCRATCHCTR_VAR,2, 0,0,0) !for scrch_Ptr - done in Hydro_prepareBuffers

  ! ***** SECOND VARIANT: FOR hy_fluxCorrectPerLevel==.FALSE. *****
  if (.NOT. hy_fluxCorrectPerLevel) then

     call Timers_start("compute fluxes")
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        level            = tileDesc%level
        blkLimits(:,:)   = tileDesc%limits
        blkLimitsGC(:,:) = tileDesc%limitsGC

        call tileDesc%getDataPtr(Uin, CENTER)

        call tileDesc%deltas(del)
        if (level==maxLev) then
           Uout => Uin             ! hy_computeFluxes will ALSO update the solution through the Uout pointer!
        else
           nullify(Uout)           ! Uout is not needed yet.
        end if
        call hy_computeFluxes(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uin, CENTER)
        nullify(Uout)
        call itor%next()
     end do
     call Timers_stop("compute fluxes")
     call Grid_releaseTileIterator(itor)

     if (hy_fluxCorrect) then
        call Grid_putFluxData(level=UNSPEC_LEVEL)
        call Timers_start("conserveFluxes")
        call Grid_conserveFluxes(ALLDIR,UNSPEC_LEVEL)
        call Timers_stop("conserveFluxes")
     end if

     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     call Timers_start("update solution")
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        level            = tileDesc%level
        if (level==maxLev) then
           call itor%next()
           CYCLE
        end if
        blkLimits(:,:)   = tileDesc%limits
        blkLimitsGC(:,:) = tileDesc%limitsGC

        call tileDesc%getDataPtr(Uout, CENTER)

        call tileDesc%deltas(del)
        Uin => Uout
        call hy_updateSolution(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uout, CENTER)
        nullify(Uin)
        call itor%next()
     end do
     call Timers_stop("update solution")
     call Grid_releaseTileIterator(itor)

!!$     call hy_memDeallocScratch(SCRATCH_CTR) !done in Hydro_freeBuffers
     RETURN                     ! DONE, return from here!
  end if


  ! ***** THIRD VARIANT: FOR hy_fluxCorrectPerLevel==.TRUE. *****

  do level= maxLev,1,-1
#ifdef DEBUG_DRIVER
     print*,' ***************   HYDRO LEVEL', level,'  **********************'
#endif
     !! if(hy_fluxCorrectPerLevel) then
     !! if(level !=maxLev) then
     !!   do a synchronization step here
     !!        if(hy_fluxCorrectPerLevel)call Grid_conserveFluxes(ALLDIR,level)
     
     call Timers_start("compute fluxes")
     call Grid_getDeltas(level, del)
     call Grid_getTileIterator(itor, LEAF, level=level, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        
        blkLimits(:,:)   = tileDesc%limits
        blkLimitsGC(:,:) = tileDesc%limitsGC
        
        call tileDesc%getDataPtr(Uin, CENTER)
        
        if (level==maxLev) then
           Uout => Uin             ! hy_computeFluxes will ALSO update the solution through the Uout pointer!
        else
           nullify(Uout)           ! Uout is not needed yet.
        end if
        call hy_computeFluxes(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uin, CENTER)
        nullify(Uout)
        call itor%next()
     end do
     call Timers_stop("compute fluxes")
     call Grid_releaseTileIterator(itor)

     if (hy_fluxCorrect .AND. (level > 1))  call Grid_putFluxData(level)

     if (level==maxLev) then
        CYCLE
     end if

     if (hy_fluxCorrect) then
        call Timers_start("conserveFluxes")
        call Grid_conserveFluxes(ALLDIR,level)
        call Timers_stop("conserveFluxes")
     end if

     call Grid_getTileIterator(itor, LEAF, level=level, tiling=.FALSE.)
     call Timers_start("update solution")
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        
        blkLimits(:,:)   = tileDesc%limits
        blkLimitsGC(:,:) = tileDesc%limitsGC
        
        call tileDesc%getDataPtr(Uout, CENTER)
 
        Uin => Uout
        call hy_updateSolution(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uout, CENTER)
        nullify(Uin)
        call itor%next()
     end do
     call Timers_stop("update solution")
     call Grid_releaseTileIterator(itor)

#ifdef DEBUG_DRIVER
     print*, 'return from Hydro/MHD timestep'  ! DEBUG
     print*,'returning from hydro myPE=',dr_globalMe
#endif


  end do
  
!!$  call hy_memDeallocScratch(SCRATCH_CTR) ! done in Hydro_freeBuffers
  
end subroutine hy_advance
