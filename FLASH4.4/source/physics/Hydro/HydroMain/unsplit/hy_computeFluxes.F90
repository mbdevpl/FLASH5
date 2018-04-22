!!****if* source/physics/Hydro/HydroMain/unsplit/hy_computeFluxes
!!
!!
!! NAME
!!
!!  hy_computeFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_computeFluxes(integer(IN) :: blockCount, 
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs physics update in a directionally unsplit fashion.
!!
!!  The blockList and blockCount arguments tell this routine on 
!!  which blocks and on how many to operate.  blockList is an 
!!  integer array of size blockCount that contains the local 
!!  block numbers of blocks on which to computeFluxes.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***


!!$#define DEBUG
!!$#define DEBUG_UHD

!!REORDER(4): scrch_Ptr, scrchFace[XYZ]Ptr, fl[xyz]

Subroutine hy_computeFluxes(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld,sweepOrder)

  use Eos_interface, ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata,   ONLY : block_metadata_t
  use hy_interface, ONLY : hy_getRiemannState,  &
                               hy_getFaceFlux,      &
                               hy_unsplitUpdate,    &
                               hy_unitConvert,      &
                               hy_energyFix,        &
                               hy_putGravity

  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_fluxCorrectPerLevel,             &
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
                         hy_fullSpecMsFluxHandling

  use Grid_interface, ONLY : Grid_getFluxPtr, Grid_releaseFluxPtr
  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real, pointer, dimension(:,:,:,:) :: Uout
  real, pointer, dimension(:,:,:,:) :: Uin

  real,dimension(MDIM),intent(IN) :: del
  integer,dimension(LOW:HIGH,MDIM),intent(INoUt) ::blkLimits,blkLimitsGC 
  integer :: loxGC,hixGC,loyGC,hiyGC,lozGC,hizGC
  type(block_metadata_t), intent(IN) :: blockDesc
  
  integer, dimension(MDIM) :: datasize

  real, pointer, dimension(:,:,:,:)   :: flx,fly,flz
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)

  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig

  integer :: updateMode ! will be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND



  call Timers_start("loop1 body")

     loxGC = blkLimitsGC(LOW,IAXIS); hixGC =blkLimitsGC(HIGH,IAXIS)
     loyGC = blkLimitsGC(LOW,JAXIS); hiyGC =blkLimitsGC(HIGH,JAXIS)
     lozGC = blkLimitsGC(LOW,KAXIS); hizGC =blkLimitsGC(HIGH,KAXIS)

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

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        call hy_unitConvert(Uin,blkLimitsGC,FWDCONVERT)
     endif

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     allocate(scrch_Ptr    (2,               loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
     allocate(scrchFaceXPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
     allocate(scrchFaceYPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
     allocate(scrchFaceZPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
!!$     endif

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        allocate(  hy_SpcR(HY_NSPEC,                                   loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC,NDIM))
        allocate(  hy_SpcL(HY_NSPEC,                                   loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC,NDIM))
        allocate(hy_SpcSig(HY_NSPEC,blkLimits(LOW,IAXIS)-2:blkLimits(HIGH,IAXIS)+2, loyGC:hiyGC, lozGC:hizGC,NDIM))
     end if
#endif

     allocate(gravX(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravY(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravZ(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
#ifdef DEBUG
     print*,'came upto this point'
#endif
     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        call hy_putGravity(blockDesc,blkLimitsGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
     endif


     if (hy_updateHydroFluxes) then
        !! ************************************************************************
        !! Calculate Riemann (interface) states
        !! Note: gravX(:,:,:) - gravity at n

#if (NSPECIES+NMASS_SCALARS) > 0
        if (hy_fullSpecMsFluxHandling) then
           hy_SpcL=0.
           hy_SpcR=0.
           hy_SpcSig=0.
        end if
#endif

#ifdef DEBUG_UHD
        print*,'_unsplit bef "call getRiemannState": associated(Uin ) is',associated(Uin )
        print*,'_unsplit bef "call getRiemannState": associated(Uout) is',associated(Uout)
        print*,'_unsplit bef "call getRiemannState": lbound(Uin ):',lbound(Uin )
        print*,'_unsplit bef "call getRiemannState": ubound(Uin ):',ubound(Uin )
        print*,'_unsplit bef "call getRiemannState": lbound(scrchFaceXPtr):',lbound(scrchFaceXPtr)
        print*,'_unsplit bef "call getRiemannState": ubound(scrchFaceXPtr):',ubound(scrchFaceXPtr)
        print*,'_unsplit bef "call getRiemannState": lbound(scrchFaceYPtr):',lbound(scrchFaceYPtr)
        print*,'_unsplit bef "call getRiemannState": ubound(scrchFaceYPtr):',ubound(scrchFaceYPtr)
#endif
        call Timers_start("RiemannState")
#ifdef DEBUG_UHD
        print*,'going into RiemannState'
#endif
        call hy_getRiemannState(blockDesc,Uin,blkLimits,blkLimitsGC(LOW,:),blkLimitsGC(HIGH,:),dt,del, &
                                    gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:),&
                                    scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                    hy_SpcR,hy_SpcL,hy_SpcSig)
#ifdef DEBUG_UHD
        print*,'returning from RiemannState'
        print*,'_unsplit Aft "call getRiemannState": associated(Uin ) is',associated(Uin )
        print*,'_unsplit Aft "call getRiemannState": associated(Uout) is',associated(Uout)
#endif
        call Timers_stop("RiemannState")
        !! DEV: DL-This note seems to be outdated and wrong for the optimized code.
        ! Note: Two different ways of handling gravity:
        ! 1. With gravity calculated potential at n+1/2, Riemann states do not include gravity
        !    source terms at this point, and will include them in hy_addGravity later
        !    to the primitive Riemann variables (not available for conservative formulation).
        ! 2. With gravity extrapolated from n-1 & n states, gravity source terms have been
        !    included to Riemann states in conservative formulation in hy_getRiemannState.

     endif !! End of if (hy_updateHydroFluxes) then
     call Grid_getFluxPtr(blockDesc,flx,fly,flz)
     allocate(  faceAreas(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))

!!$     call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR) 

     !! ************************************************************************
     !! Calculate high order Godunov fluxes
     !! Initialize arrays with zero
     flx = 0.
     fly = 0.
     flz = 0.
     call Timers_start("getFaceFlux")
#ifdef DEBUG_UHD
     print*,'getting face flux'
#endif
     call hy_getFaceFlux(blockDesc,blkLimits,blkLimitsGC,datasize,del,&
                             flx,fly,flz,&
                             scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,hy_SpcR,hy_SpcL)
#ifdef DEBUG_UHD
     print*,'got face flux'
     print*,'_unsplit Aft "call getFaceFlux": associated(Uin ) is',associated(Uin )
     print*,'_unsplit Aft "call getFaceFlux": associated(Uout) is',associated(Uout)
#endif
     call Timers_stop("getFaceFlux")
     !! ************************************************************************
     !! Unsplit update for conservative variables from n to n+1 time step

     
     call Grid_releaseFluxPtr(blockDesc,flx,fly,flz)
     
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
     
     deallocate(scrchFaceXPtr)
     deallocate(scrchFaceYPtr)
     deallocate(scrchFaceZPtr)
     
#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        deallocate(hy_SpcR)
        deallocate(hy_SpcL)
        deallocate(hy_SpcSig)
     end if
#endif

  call Timers_stop("loop1 body")


End Subroutine hy_computeFluxes
