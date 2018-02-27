#define DEBUG_GRID_GCMASK

#include "constants.h"


subroutine Hydro_advanceAll(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_fillGuardCells
  use Grid_interface,      ONLY : Grid_getMaxRefinement
  use Grid_interface,      ONLY : Grid_copyF4DataToMultiFabs
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_prepareBuffers, Hydro_freeBuffers
  use Hydro_interface,     ONLY : Hydro_shockDetectLoop, Hydro_computeFluxLoop, Hydro_doLoop4
  use Hydro_data, ONLY : hy_gcMaskSize,       &
                         hy_gcMask

#include "Flash.h"
#ifdef FLASH_GRID_AMREXTRANSITION
  use gr_amrextInterface,  ONLY : gr_amrextBuildMultiFabsFromF4Grid
  use gr_amrextData
#endif

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  integer:: maxLev

  call Timers_start("Hydro")

#ifdef FLASH_GRID_AMREXTRANSITION
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif

  call Hydro_prepareBuffers()

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[Hydro_advanceAll]', 'gcNeed')
  end if
#endif
  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
       maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)

#ifdef FLASH_GRID_AMREXTRANSITION
  call gr_amrextBuildMultiFabsFromF4Grid(gr_amrextUnkMFs, maxLev, LEAF)
#endif
  call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)

  call Hydro_shockDetectLoop()

  call Hydro_computeFluxLoop(simTime, dt, dtOld)

  call Hydro_doLoop4()

  call Hydro_freeBuffers()

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  call Timers_stop("Hydro")

end subroutine Hydro_advanceAll
