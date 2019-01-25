#include "Flash.h"
#include "constants.h"

subroutine hy_gravityStep(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getTileIterator, &
                                  Grid_releaseTileIterator
  use Timers_interface,    ONLY : Timers_start, Timers_stop

  use hy_interface,        ONLY : hy_gravityStepBlk
  use flash_iterator,      ONLY : flash_iterator_t
  use flash_tile,          ONLY : flash_tile_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

#ifdef DEBUG_DRIVER
  print*,' ***************   HYDRO LEVEL', level,'  **********************'
#endif

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  call Timers_stop("loop5")
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     blkLimits(:,:)   = tileDesc%limits
     blkLimitsGC(:,:) = tileDesc%limitsGC
     
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
