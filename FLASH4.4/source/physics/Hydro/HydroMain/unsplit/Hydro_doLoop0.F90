#include "constants.h"

subroutine Hydro_doLoop0
  use Grid_interface, ONLY : Grid_getDeltas, &
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_prepareNewGravityAccel,&
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_shockDetect
  use block_iterator, ONLY : block_iterator_t, destroy_iterator
  use block_metadata, ONLY : block_metadata_t
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
  implicit none

#include "UHD.h"

  real, dimension(MDIM) :: del

  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:),pointer :: Uin
  real,dimension(:,:,:,:), pointer :: Uout


!!$  do i=1,blockCount          !LOOP 0
!!$     blockID = blockList(i)

  itor = block_iterator_t(LEAF)

  do while(itor%is_valid())
     call itor%blkMetaData(blockDesc)

     blkLimits(:,:)   = blockDesc%localLimits
     blkLimitsGC(:,:) = blockDesc%localLimitsGC

     call Grid_getBlkPtr(blockDesc, Uout,localFlag=.TRUE.)
     Uin => Uout

     !! Detect shocks
     if (hy_shockDetectOn) then
        call Grid_getDeltas(blockDesc%level,del)
        call hy_uhd_shockDetect(Uin,blkLimitsGC,Uout,blkLimits,del)
     end if

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        call hy_uhd_unitConvert(Uin,blkLimitsGC,FWDCONVERT)
     endif
     
#if defined(GPRO_VAR)||defined(VOLX_VAR)||defined(VOLY_VAR)||defined(VOLZ_VAR)||defined(CFL_VAR)
     if (hy_updateHydroFluxes) then
#ifdef GPRO_VAR
        ! A tagging variable for Gaussian Process (GP) method.
        Uin(GPRO_VAR,:,:,:) = 0.
#endif
        !! -----------------------------------------------------------------------!
        !! Save old velocities ---------------------------------------------------!
        !! -----------------------------------------------------------------------!
#ifdef VOLX_VAR
        Uin(VOLX_VAR,:,:,:) = Uin(VELX_VAR,:,:,:)
#endif
#ifdef VOLY_VAR
        Uin(VOLY_VAR,:,:,:) = Uin(VELY_VAR,:,:,:)
#endif
#ifdef VOLZ_VAR
        Uin(VOLZ_VAR,:,:,:) = Uin(VELZ_VAR,:,:,:)
#endif
#ifdef CFL_VAR
        where (1.2*Uin(CFL_VAR,:,:,:) < hy_cfl_original)
           !! Slow recover (of factor of 1.2) to the original CFL once it gets to
           !! reduced to a smaller one in the presence of strong shocks.
           !! This variable CFL takes place in the following three cases using:
           !! (1) use_hybridOrder = .true.,
           !! (2) use_hybridOrder = .true., or
           !! (3) BDRY_VAR is defined and used for stationary objects.
           Uin(CFL_VAR,:,:,:) = 1.2*U(CFL_VAR,:,:,:)
        elsewhere
           Uin(CFL_VAR,:,:,:) = hy_cfl_original
        end where
#endif
     end if
#endif
     call Grid_releaseBlkPtr(blockDesc, Uout)

     call itor%next()
  end do
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
  call destroy_iterator(itor)
#endif
end subroutine Hydro_doLoop0
