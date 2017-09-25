#define DEBUG_GRID_GCMASK

#include "constants.h"


subroutine Hydro_advanceAll(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_fillGuardCells
  use Grid_interface,      ONLY : Grid_copyF4DataToMultiFabs
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_prepareBuffers, Hydro_freeBuffers
  use Hydro_interface,     ONLY : Hydro_doLoop0, Hydro_doLoop1, Hydro_doLoop4
  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_units,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,    &
                         hy_eosModeGc,        &
                         hy_eosModeAfter,     &
                         hy_updateHydroFluxes,&
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_dtmin,            &
                         hy_simTime,          &
                         hy_simGeneration,    &
                         hy_shockDetectOn,    &
                         hy_doUnsplitLoop0
  use Hydro_data,       ONLY : hy_useHydro, hy_gpotAlreadyUpToDate

#include "Flash.h"
#ifdef FLASH_GRID_AMREXTRANSITION
  use gr_amrextInterface,  ONLY : gr_amrextBuildMultiFabsFromF4Grid
!!$  use gr_amrextData
#endif

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  logical :: gcMask(hy_gcMaskSize)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  hy_gpotAlreadyUpToDate = .FALSE. ! reset this flag, may be set .TRUE. below if warranted.

  if (.not. hy_useHydro) return

  call Timers_start("Hydro")

  call Timers_start("Head")

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The unsplit Hydro solver only works with PARAMESH 3 or 4!")
#endif

#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
#endif

  call Hydro_prepareBuffers()

#ifdef FLASH_GRID_AMREXTRANSITION
     call gr_amrextBuildMultiFabsFromF4Grid(gr_amrextUnkMFs, maxLev, LEAF)
#endif
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)


  !! ***************************************************************************
  !! Shock detection before hydro (etc.)                                       *
  !! ***************************************************************************
  !! Call shock detect algorithm to determine tagging shocks before hydro begins;
  !! other preparations of UNK data if necessary.
  if (hy_shockDetectOn) then
     
     !! Call guardcell filling to properly detect shocks
     gcMask = .false.
     gcMask(DENS_VAR) = .true.
     gcMask(PRES_VAR) = .true.
     gcMask(GAMC_VAR) = .true.
     gcMask(VELX_VAR:VELZ_VAR) = .true.
#ifdef CFL_VAR
     gcMask(CFL_VAR)  = .true.
#endif
#if NSPECIES > 1
     gcMask(SPECIES_BEGIN:SPECIES_END) = .true.
#endif
     
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Detect]')
     end if
#endif
     
     call Grid_fillGuardCells(CENTER,ALLDIR & !) ! DEV: NONONO!
     ,doEos=.false.,&
          maskSize=NUNK_VARS, mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
     
     hy_cfl = hy_cfl_original
  end if


  if (hy_doUnsplitLoop0) then
     call Hydro_doLoop0()
  endif


  !! ***************************************************************************
  !! Call guardcell filling with Eos before hydro                              *
  !! ***************************************************************************
#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcNeed')
  end if
#endif

  ! DEV: NONONO!
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
  !! Loop over the blocks
  call Hydro_doLoop1(simTime, dt, dtOld)
!!$  call IO_writeCheckpoint()
!!$  stop

  call Hydro_doLoop4()



  call Hydro_freeBuffers()



  call Driver_getSimTime(hy_simTime, hy_simGeneration)

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  call Timers_stop("Hydro")

end subroutine Hydro_advanceAll
