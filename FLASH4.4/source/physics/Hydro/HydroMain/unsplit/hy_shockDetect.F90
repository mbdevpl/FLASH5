#include "constants.h"

subroutine hy_shockDetect
!  use Grid_interface, ONLY : Grid_getDeltas, &
!                             Grid_getBlkPtr, Grid_releaseBlkPtr, &
!                             Grid_getLeafIterator, Grid_releaseLeafIterator
  use hy_interface,   ONLY : hy_shockDetectBlk
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_units,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeGc,        &
                         hy_eosModeAfter,     &
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_numXN,            &
                         hy_fullRiemannStateArrays,    &
                         hy_fullSpecMsFluxHandling,   &
                         hy_dtmin,            &
                         hy_simTime,          &
                         hy_simGeneration,    &
                         hy_shockDetectOn
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "UHD.h"

!  real, dimension(MDIM) :: del
!
!  type(Grid_iterator_t) :: itor
!  type(Grid_tile_t)     :: tileDesc
!
!  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
!  real, dimension(:,:,:,:),pointer :: Uin
!  real,dimension(:,:,:,:), pointer :: Uout

  call Driver_abortFlash("[hy_shockDetect] Implement for use with tiling")

!  call Grid_getLeafIterator(itor)
!  do while(itor%is_valid())
!     call itor%blkMetaData(blockDesc)
!     
!     blkLimits(:,:)   = blockDesc%localLimits
!     blkLimitsGC(:,:) = blockDesc%localLimitsGC
!     
!     call Grid_getBlkPtr(blockDesc, Uout,localFlag=.TRUE.)
!     Uin => Uout
!     
!     !! Detect shocks
!     if (hy_shockDetectOn) then
!        call Grid_getDeltas(blockDesc%level,del)
!        call hy_shockDetectBlk(Uin,blkLimitsGC,Uout,blkLimits,del)
!     end if
!     
!     
!     call Grid_releaseBlkPtr(blockDesc, Uout)
!     
!     call itor%next()
!  end do
!  call Grid_releaseLeafIterator(itor)
  
end subroutine hy_shockDetect
