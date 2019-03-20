#include "constants.h"

subroutine hy_shockDetect
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Hydro_data, ONLY : hy_shockDetectOn
  use hy_interface,   ONLY : hy_shockDetectBlk

  implicit none

#include "UHD.h"

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  integer, dimension(LOW:HIGH,MDIM) :: limits,grownLimits

  real, dimension(MDIM) :: del

  real, dimension(:,:,:,:),pointer :: Uin
  real,dimension(:,:,:,:), pointer :: Uout

  !! DEV: All executabe statements that follow could be skipped
  !! completely if the hy_shockDetectOn flag is TRUE; this is not
  !! done here as a reminder that the corresponding loop in FLASH4
  !! code [LOOP 0 in hy_uhd_unsplit.F90] was doing more than just
  !! shock detection.

  call Grid_getTileIterator(itor,LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     limits(:,:)   = tileDesc%limits
     grownLimits(:,:) = tileDesc%grownLimits

     call tileDesc % getDataPtr(Uout,CENTER)
     Uin => Uout

     !! Detect shocks
     if (hy_shockDetectOn) then
        call tileDesc % deltas(del)
        call hy_shockDetectBlk(Uin,lbound(Uin),grownLimits,Uout,lbound(Uout),limits,del)
     end if


     call tileDesc % releaseDataPtr(Uout,CENTER)

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  
end subroutine hy_shockDetect
