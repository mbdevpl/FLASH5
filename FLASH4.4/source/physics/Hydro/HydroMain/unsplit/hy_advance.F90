#include "Flash.h"
#include "constants.h"

subroutine hy_advance(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getDeltas,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_getBlkIterator, Grid_releaseBlkIterator,&
                                  Grid_getMaxRefinement
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use hy_interface,        ONLY : hy_advanceBlk
  use block_iterator, ONLY : block_iterator_t
  use block_metadata, ONLY : block_metadata_t
  use Hydro_data, ONLY : hy_fluxCorrectPerLevel

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
        !! if(hy_fluxCorrectPerLevel) then
        !! if(level !=maxLev) then
        !!   do a synchronization step here
        if(hy_fluxCorrectPerLevel)call Grid_conserveFluxes(ALLDIR,level)
        call Grid_getBlkIterator(itor, LEAF, level=level)
        call Timers_stop("loop1")
        do while(itor%is_valid())
           call itor%blkMetaData(blockDesc)

           blkLimits(:,:)   = blockDesc%limits
           blkLimitsGC(:,:) = blockDesc%limitsGC
           
           call Grid_getBlkPtr(blockDesc, Uout)
!!$           abx = amrex_box(bx%lo, bx%hi, bx%nodal)
!!$           call amrex_print(abx)
!!$           tbx = abx

           call Grid_getDeltas(level,del)
           Uin => Uout
           call hy_advanceBlk(blockDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld,  sweepDummy)
           call Grid_releaseBlkPtr(blockDesc, Uout)
           nullify(Uout)
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


end subroutine hy_advance
