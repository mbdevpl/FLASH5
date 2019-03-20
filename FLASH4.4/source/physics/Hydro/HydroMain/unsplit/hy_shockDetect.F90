#include "constants.h"

subroutine hy_shockDetect
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use hy_interface,   ONLY : hy_shockDetectBlk
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Hydro_data, ONLY : hy_shockDetectOn
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "UHD.h"

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
!
  integer, dimension(LOW:HIGH,MDIM) :: limits,grownLimits

  real, dimension(MDIM) :: del
!
  real, dimension(:,:,:,:),pointer :: Uin
  real,dimension(:,:,:,:), pointer :: Uout

  call Grid_getTileIterator(itor,LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
!     
     limits(:,:)   = tileDesc%limits
     grownLimits(:,:) = tileDesc%grownLimits
!     
     call tileDesc % getDataPtr(Uout,CENTER,localFlag=.FALSE.)
     Uin => Uout
!     
     !! Detect shocks
     if (hy_shockDetectOn) then
        call tileDesc % deltas(del)
        call hy_shockDetectBlk(Uin,lbound(Uin),grownLimits,Uout,lbound(Uout),limits,del)
!        call hy_shockDetectBlk(Uin,blkLimitsGC,Uout,blkLimits,del)
     end if
!     
!     
     call tileDesc % releaseDataPtr(Uout,CENTER)
!     
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  
end subroutine hy_shockDetect
