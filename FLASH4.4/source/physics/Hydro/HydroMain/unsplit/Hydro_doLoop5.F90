#include "Flash.h"
#include "constants.h"

subroutine Hydro_doLoop5(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getDeltas,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_getMaxRefinement
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Hydro_interface,     ONLY : Hydro_loop5Body
  use block_iterator, ONLY : block_iterator_t
  use block_metadata, ONLY : block_metadata_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  integer:: level, maxLev

  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.

  do level=1,maxLev
#ifdef DEBUG_DRIVER
        print*,' ***************   HYDRO LEVEL', level,'  **********************'
#endif

        itor = block_iterator_t(LEAF, level=level)
        call Timers_stop("loop5")
        do while(itor%is_valid())
           call itor%blkMetaData(blockDesc)

           blkLimits(:,:)   = blockDesc%limits
           blkLimitsGC(:,:) = blockDesc%limitsGC
           
           call Grid_getBlkPtr(blockDesc, Uout)

           call Grid_getDeltas(level,del)
           Uin => Uout
           call Hydro_loop5Body(blockDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld)
           call Grid_releaseBlkPtr(blockDesc, Uout)
           nullify(Uout)
!!$           call IO_writecheckpoint;stop
           call itor%next()
        end do
        call Timers_stop("loop5")
        
        
     end do


end subroutine Hydro_doLoop5
