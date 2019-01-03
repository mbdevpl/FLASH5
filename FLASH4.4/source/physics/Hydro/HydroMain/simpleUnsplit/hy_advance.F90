!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/hy_advance
!!
!! NAME
!!  hy_advance
!!
!! DESCRIPTION
!!  Refer to stub for detailed documentation.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine hy_advance(simTime, dt, dtOld)
  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop
  use hy_interface,     ONLY : hy_advanceBlk
  use flash_iterator,   ONLY : flash_iterator_t
  use flash_tile,       ONLY : flash_tile_t

  implicit none

  real, intent(IN) :: simTime
  real, intent(IN) :: dt
  real, intent(IN) :: dtOld

  real, pointer :: Uout(:,:,:,:)
  real, pointer :: Uin(:,:,:,:)

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  call Timers_stop("loop1")
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%getDataPtr(Uout, CENTER)
     Uin => Uout

     call hy_advanceBlk(tileDesc, Uin, Uout, simTime, dt, dtOld, SWEEP_ALL)

     call tileDesc%releaseDataPtr(Uout, CENTER)
 
     call itor%next()
  end do
  call Timers_stop("loop1")
  call Grid_releaseTileIterator(itor)

#ifdef DEBUG_DRIVER
  print*, 'return from Hydro/MHD timestep'  ! DEBUG
  print*,'returning from hydro myPE=',dr_globalMe
#endif

end subroutine hy_advance

