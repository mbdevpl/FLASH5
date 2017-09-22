#include "constants.h"


subroutine Hydro_advanceAll(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getMaxRefinement
  use Grid_interface,      ONLY : Grid_copyF4DataToMultiFabs
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_prepareBuffers, Hydro_freeBuffers
  use Hydro_interface,     ONLY : Hydro_doLoop0, Hydro_doLoop1, Hydro_doLoop4

#include "Flash.h"
#ifdef FLASH_GRID_AMREXTRANSITION
  use gr_amrextInterface,  ONLY : gr_amrextBuildMultiFabsFromF4Grid
  use gr_amrextData
#endif

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer:: maxLev

  call Timers_start("Hydro")

#ifdef FLASH_GRID_AMREXTRANSITION
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif

  call Hydro_prepareBuffers()

#ifdef FLASH_GRID_AMREXTRANSITION
     call gr_amrextBuildMultiFabsFromF4Grid(gr_amrextUnkMFs, maxLev, LEAF)
#endif
     call Grid_copyF4DataToMultiFabs(CENTER, nodetype=LEAF)

     call Hydro_doLoop0()

     call Hydro_doLoop1(simTime, dt, dtOld)

     call Hydro_doLoop4()

     call Hydro_freeBuffers()

  call Timers_stop("Hydro")

end subroutine Hydro_advanceAll
