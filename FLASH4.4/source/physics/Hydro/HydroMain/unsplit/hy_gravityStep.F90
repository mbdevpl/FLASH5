#include "Flash.h"
#include "constants.h"

subroutine hy_gravityStep(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getTileIterator, &
                                  Grid_releaseTileIterator
  use Timers_interface,    ONLY : Timers_start, Timers_stop

  use hy_interface,        ONLY : hy_gravityStepBlk
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Grid_tile,           ONLY : Grid_tile_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  nullify(Uin)
  nullify(Uout)

#ifdef DEBUG_DRIVER
  print*,' ***************   HYDRO LEVEL  **********************'
#endif

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  call Timers_stop("loop5")
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     blkLimits(:,:)   = tileDesc%limits
     blkLimitsGC(:,:) = tileDesc%blkLimitsGC

     call tileDesc%getDataPtr(Uout, CENTER)
     call tileDesc%deltas(del)
     Uin => Uout
     call hy_gravityStepBlk(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld)
     call tileDesc%releaseDataPtr(Uout, CENTER)
     nullify(Uout)
!!      call IO_writecheckpoint;stop
     call itor%next()
  end do
  call Timers_stop("loop5")
  call Grid_releaseTileIterator(itor)
end subroutine hy_gravityStep
