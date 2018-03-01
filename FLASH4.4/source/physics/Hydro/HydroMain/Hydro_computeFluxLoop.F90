#include "Flash.h"
#include "constants.h"

subroutine Hydro_computeFluxLoop(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getDeltas,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_getBlkIterator, Grid_releaseBlkIterator,&
                                  Grid_getMaxRefinement, Grid_conserveFluxes
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use hy_interface,        ONLY : hy_computeFluxes
  use block_iterator, ONLY : block_iterator_t
  use block_metadata, ONLY : block_metadata_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  integer,save :: sweepDummy = SWEEP_ALL

  integer:: level, maxLev

  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.

  do level=1,maxLev
#ifdef DEBUG_DRIVER
        print*,' ***************   HYDRO LEVEL', level,'  **********************'
#endif

        call Grid_getBlkIterator(itor, LEAF, level=level)
        call Timers_stop("loop1")
        do while(itor%is_valid())
           call itor%blkMetaData(blockDesc)
           
           blkLimits(:,:)   = blockDesc%limits
           blkLimitsGC(:,:) = blockDesc%limitsGC
           
!!$ DEV-AD the burden of having separate uin and uout managed here, perhaps an additional argument needed. 
!!$ Perhaps Uout can be made to be just blocksize, the question is who should have the intelligence to 
!!$ to return the right pointer.
           call Grid_getBlkPtr(blockDesc, Uout)
!!$           call Grid_getBlkPtr(blockDesc, Uin)
           
!!$ DEV-AD Here we get the storage for the faces where computed fluxes will be stored.
!!$ DEV-AD call Grid_getFaceBlkPtr(blockDesc,flxx, flxy, flxz)

!!$           abx = amrex_box(bx%lo, bx%hi, bx%nodal)
!!$           call amrex_print(abx)
!!$           tbx = abx
           
           call Grid_getDeltas(level,del)
           Uin => Uout
           call hy_advance(blockDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)

!!$           call Grid_conserveFluxes(ALLDIR,level)
!!$           call hy_advance(blockDesc,blkLimitsGC,blkLimits, uin, Uout, flxx, flxy, flxz,&
!!$                               del,simTime, dt, dtOld,  sweepDummy)

           call Grid_releaseBlkPtr(blockDesc, Uout)
           nullify(Uout)
!!$           call Grid_releaseBlkPtr(blockDesc, Uin)


!!$           call IO_writecheckpoint;stop
           call itor%next()
        end do
        call Timers_stop("loop1")
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
        call Grid_releaseBlkIterator(itor)
#endif
#ifdef DEBUG_DRIVER
        print*, 'return from Hydro/MHD timestep'  ! DEBUG
        print*,'returning from hydro myPE=',dr_globalMe
#endif
        
        
     end do


end subroutine Hydro_computeFluxLoop
