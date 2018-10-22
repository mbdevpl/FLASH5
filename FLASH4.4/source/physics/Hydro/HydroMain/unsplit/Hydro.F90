#define DEBUG_GRID_GCMASK

#include "constants.h"


subroutine Hydro(simTime, dt, dtOld, sweeporder)

  use Grid_interface,      ONLY : Grid_fillGuardCells
  use Grid_interface,      ONLY : Grid_getMaxRefinement
  use Grid_interface,      ONLY : Grid_copyF4DataToMultiFabs
  use Driver_interface,    ONLY : Driver_getSimTime
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_prepareBuffers, Hydro_freeBuffers
  use Hydro_data, ONLY : hy_fluxCorrect,       &
                         hy_fluxCorrectPerLevel, &
                         hy_gref,              &
                         hy_useGravity,        &
                         hy_units,             &
                         hy_gcMaskSize,        &
                         hy_gcMask,hy_gcMaskSD,&
                         hy_eosModeGc,         &
                         hy_eosModeAfter,      &
                         hy_updateHydroFluxes, &
                         hy_cfl,               &
                         hy_cfl_original,      &
                         hy_dtmin,             &
                         hy_simTime,           &
                         hy_simGeneration,     &
                         hy_shockDetectOn

  use Hydro_data,       ONLY : hy_useHydro, hy_gpotAlreadyUpToDate
  use hy_interface, ONLY : hy_gravityStep

#include "Flash.h"
#ifdef FLASH_GRID_AMREXTRANSITION
  use gr_amrextInterface,  ONLY : gr_amrextBuildMultiFabsFromF4Grid
#endif

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld
  integer, optional, intent(IN):: sweeporder

!!  logical :: gcMask(hy_gcMaskSize)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  integer:: maxLev

  hy_gpotAlreadyUpToDate = .FALSE. ! reset this flag, may be set .TRUE. below if warranted.

  if (.not. hy_useHydro) return

  call Timers_start("Hydro")

  call Timers_start("Head")

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The unsplit Hydro solver only works with PARAMESH 3 or 4!")
#endif

#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
  maxLev = 1
#else
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
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


#ifdef FLASH_GRID_AMREXTRANSITION
     call gr_amrextBuildMultiFabsFromF4Grid(CENTER, maxLev, LEAF)
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)
#endif
     call hy_shockDetect()
#ifdef FLASH_GRID_AMREXTRANSITION
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF,reverse=.TRUE.)
     call gr_amrextBuildMultiFabsFromF4Grid(CENTER, maxLev, ACTIVE_BLKS)
#endif
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
#ifdef FLASH_GRID_AMREXTRANSITION
  call gr_amrextBuildMultiFabsFromF4Grid(CENTER, maxLev, LEAF)
  call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)
#endif


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
  call hy_advance(simTime, dt, dtOld)

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
#ifdef FLASH_GRID_AMREXTRANSITION
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF,reverse=.TRUE.)
     call gr_amrextBuildMultiFabsFromF4Grid(CENTER, maxLev, ACTIVE_BLKS)
#endif
     ! The following call invokes Gravity_potential and related stuff,
     ! to prepare for retrieving updated accelerations below.
     call hy_prepareNewGravityAccel(gcMaskLogged)
#ifdef FLASH_GRID_AMREXTRANSITION
     call gr_amrextBuildMultiFabsFromF4Grid(CENTER, maxLev, LEAF)
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)
#endif
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
