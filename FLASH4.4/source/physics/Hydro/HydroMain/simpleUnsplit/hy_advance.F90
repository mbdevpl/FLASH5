#include "Flash.h"
#include "constants.h"

subroutine hy_advance(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getDeltas,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_getLeafIterator, Grid_releaseLeafIterator
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use hy_interface,    ONLY : hy_advanceBlk
  use leaf_iterator,       ONLY : leaf_iterator_t
  use block_metadata,      ONLY : block_metadata_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  integer,save :: sweepDummy = SWEEP_ALL

  type(leaf_iterator_t)  :: itor
  type(block_metadata_t) :: blockDesc

  call Grid_getLeafIterator(itor, tiling=.FALSE.)
  call Timers_stop("loop1")
  do while(itor%is_valid())
     call itor%blkMetaData(blockDesc)

     blkLimits(:,:)   = blockDesc%limits(:,:)
     blkLimitsGC(:,:) = blockDesc%limitsGC(:,:)
     
     call Grid_getBlkPtr(blockDesc, Uout, localFlag=.FALSE.)

     call Grid_getDeltas(blockDesc%level,del)
     Uin => Uout
     call hy_advanceBlk(blockDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
     call Grid_releaseBlkPtr(blockDesc, Uout)
 
     call itor%next()
  end do
  call Timers_stop("loop1")
  call Grid_releaseLeafIterator(itor)

#ifdef DEBUG_DRIVER
  print*, 'return from Hydro/MHD timestep'  ! DEBUG
  print*,'returning from hydro myPE=',dr_globalMe
#endif

end subroutine hy_advance

