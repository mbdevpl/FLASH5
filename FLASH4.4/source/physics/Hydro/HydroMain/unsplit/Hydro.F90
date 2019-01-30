#define DEBUG_GRID_GCMASK

#include "Flash.h"
#include "constants.h"

subroutine Hydro(simTime, dt, dtOld, sweeporder)

  use Grid_interface,      ONLY : Grid_fillGuardCells, &
                                  Grid_getMaxRefinement, &
                                  Grid_getDeltas, &
                                  Grid_zeroFluxData, &
                                  Grid_getTileIterator, &
                                  Grid_releaseTileIterator, &
                                  Grid_conserveFluxes, &
                                  Grid_putFluxData
  use Grid_data,           ONLY : gr_enableTiling
  use Driver_interface,    ONLY : Driver_getSimTime
  use Logfile_interface,   ONLY : Logfile_stampVarMask
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_prepareBuffers, Hydro_freeBuffers
  use Hydro_data,          ONLY : hy_fluxCorrect,       &
                                  hy_fluxCorrectPerLevel, &
                                  hy_useGravity,        &
                                  hy_gcMaskSize,        &
                                  hy_gcMask,hy_gcMaskSD,&
                                  hy_eosModeGc,         &
                                  hy_updateHydroFluxes, &
                                  hy_cfl,               &
                                  hy_cfl_original,      &
                                  hy_dtmin,             &
                                  hy_simTime,           &
                                  hy_simGeneration,     &
                                  hy_shockDetectOn,     &
                                  hy_useHydro,          &
                                  hy_gpotAlreadyUpToDate
  use hy_interface,        ONLY : hy_gravityStep, &
                                  hy_computeFluxes, &
                                  hy_updateSolution
  use flash_iterator,      ONLY : flash_iterator_t
  use flash_tile,          ONLY : flash_tile_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld
  integer, optional, intent(IN):: sweeporder

  integer, save :: sweepDummy = SWEEP_ALL
  
!!  logical :: gcMask(hy_gcMaskSize)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  integer:: level, maxLev
  
  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc
  
  real, pointer :: Uin(:,:,:,:)
  real, pointer :: Uout(:,:,:,:)
  real :: del(1:MDIM)

  logical :: useTiling

  nullify(Uin)
  nullify(Uout)

  hy_gpotAlreadyUpToDate = .FALSE. ! reset this flag, may be set .TRUE. below if warranted.

  if (.not. hy_useHydro) return

  call Timers_start("Hydro")

  call Timers_start("Head")

#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
  maxLev = 1
#else
  ! mode=1 means lrefine_max, which does not change during sim.
  call Grid_getMaxRefinement(maxLev, mode=1) 
#endif

  call Hydro_prepareBuffers()

  !! ***************************************************************************
  !! Shock detection before hydro (etc.)                                       *
  !! ***************************************************************************
  !! Call shock detect algorithm to determine tagging shocks before hydro begins;
  !! other preparations of UNK data if necessary.
  if (hy_shockDetectOn) then
     
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(hy_gcMaskSD, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Detect]')
     end if
#endif
     
     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,&
          maskSize=NUNK_VARS, mask=hy_gcMaskSD,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
     
     hy_cfl = hy_cfl_original


     call hy_shockDetect()
  endif

  !! ***************************************************************************
  !! Call guardcell filling with Eos before hydro                              *
  !! ***************************************************************************
#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcNeed')
  end if
#endif

  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,eosMode=hy_eosModeGc,&
       maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)

  call Timers_stop("Head")


  !! Retain the original cfl that may have been changed in some leaf blocks.
  if (hy_updateHydroFluxes) then
     if (1.2*hy_cfl < hy_cfl_original) then
        !! Slow recover (of factor of 1.2) to the original CFL once it gets to
        !! reduced to a smaller one in the presence of strong shocks.
        !! This variable CFL takes place in the following three cases using:
        !! (1) use_hybridOrder = .true.,
        !! (2) use_hybridOrder = .true., or
        !! (3) BDRY_VAR is defined and used for stationary objects.
        hy_cfl = 1.2*hy_cfl
     else
        hy_cfl = hy_cfl_original
     endif
     hy_dtmin = huge(1.0)
  endif

  !! ***************************************************************************
  !! First part of advancement                                                 *
  !! ***************************************************************************
  call Grid_zeroFluxData

  if (.NOT. hy_fluxCorrect) then
     ! ***** FIRST VARIANT: OPTIMIZED (somewhat) FOR hy_fluxCorrect==.FALSE. *****
     
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     call Timers_start("advance")
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        call tileDesc%getDataPtr(Uin, CENTER)
        call tileDesc%deltas(del)
        Uout => Uin             ! hy_computeFluxes will ALSO update the solution through the Uout pointer!
        call hy_computeFluxes(tileDesc, Uin, Uout, del, simTime, dt, dtOld, sweepDummy)
        call tileDesc%releaseDataPtr(Uin, CENTER)
        nullify(Uout)

        call itor%next()
     end do
     call Timers_stop("advance")
     call Grid_releaseTileIterator(itor)
  else if (.NOT. hy_fluxCorrectPerLevel) then
     ! ***** SECOND VARIANT: FOR hy_fluxCorrectPerLevel==.FALSE. *****

     call Timers_start("compute fluxes")
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        level = tileDesc%level
        call tileDesc%deltas(del)
        call tileDesc%getDataPtr(Uin, CENTER)
        if (level==maxLev) then
           Uout => Uin             ! hy_computeFluxes will ALSO update the solution through the Uout pointer!
        else
           nullify(Uout)           ! Uout is not needed yet.
        end if
        call hy_computeFluxes(tileDesc, Uin, Uout, del, simTime, dt, dtOld, sweepDummy)
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

        level = tileDesc%level
        if (level==maxLev) then
           call itor%next()
           CYCLE
        end if

        call tileDesc%getDataPtr(Uout, CENTER)

        call tileDesc%deltas(del)
        Uin => Uout
        call hy_updateSolution(tileDesc,Uin, Uout, del,simTime, dt, dtOld,  sweepDummy)
        call tileDesc%releaseDataPtr(Uout, CENTER)
        nullify(Uin)

        call itor%next()
     end do
     call Timers_stop("update solution")
     call Grid_releaseTileIterator(itor)
  else
     ! ***** THIRD VARIANT: FOR hy_fluxCorrectPerLevel==.TRUE. *****
     useTiling = gr_enableTiling
     do level= maxLev,1,-1
#ifdef DEBUG_DRIVER
        print*,' ***************   HYDRO LEVEL', level,'  **********************'
#endif

        call Timers_start("compute fluxes")
        call Grid_getDeltas(level, del)
        call Grid_getTileIterator(itor, LEAF, level=level, tiling=useTiling)
        do while(itor%isValid())
           call itor%currentTile(tileDesc)

           call tileDesc%getDataPtr(Uin, CENTER)
           if ((level==maxLev) .AND. (.NOT. useTiling)) then
              ! hy_computeFluxes will ALSO update the solution through the Uout pointer! 
              ! This is not compatible with tiling/stencil-based computations as
              ! the values computed in the interior of some tiles will use values
              ! already computed for this time step as opposed to from the last step
              Uout => Uin
           else
               ! Otherwise, no need to store solutions yet
              nullify(Uout)
           end if
           call hy_computeFluxes(tileDesc, Uin, Uout, &
                                 del, simTime, dt, dtOld, sweepDummy)
           call tileDesc%releaseDataPtr(Uin, CENTER)
           nullify(Uout)

           call itor%next()
        end do
        call Timers_stop("compute fluxes")
        call Grid_releaseTileIterator(itor)

        if (hy_fluxCorrect .AND. (level > 1))  call Grid_putFluxData(level)

        if ((level==maxLev) .AND. (.NOT. useTiling)) then
           ! We already have the updated solution in this special, optimized case
           ! and there is no need to do flux correction.
           CYCLE
        end if

        if (hy_fluxCorrect) then
           call Timers_start("conserveFluxes")
           call Grid_conserveFluxes(ALLDIR,level)
           call Timers_stop("conserveFluxes")
        end if

        call Grid_getTileIterator(itor, LEAF, level=level, tiling=.TRUE.)
        call Timers_start("update solution")
        do while(itor%isValid())
           call itor%currentTile(tileDesc)

           call tileDesc%getDataPtr(Uout, CENTER)
           Uin => Uout
           call hy_updateSolution(tileDesc, Uin, Uout, &
                                  del, simTime, dt, dtOld, sweepDummy)
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
  end if
  
  if(.not.hy_fluxCorrectPerLevel)then
     call hy_updateBoundaries()
  end if

  call Hydro_freeBuffers()



#ifdef GRAVITY /* Perform this only when gravity is used */
  !! ***************************************************************************
  !! Fourth part of advancement to compute gravity at n+1 state                *
  !! ***************************************************************************

#ifdef GPOT_VAR
  if (hy_useGravity) then
     ! The following call invokes Gravity_potential and related stuff,
     ! to prepare for retrieving updated accelerations below.
     call hy_prepareNewGravityAccel(gcMaskLogged)
  endif
#endif

  call hy_gravityStep(simTime, dt, dtOld)

#endif /* End of n+1 gravity coupling */



  call Driver_getSimTime(hy_simTime, hy_simGeneration)

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  call Timers_stop("Hydro")

end subroutine Hydro
