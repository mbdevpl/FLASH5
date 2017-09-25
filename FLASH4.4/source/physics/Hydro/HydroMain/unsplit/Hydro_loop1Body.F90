!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_loop1Body
!!
!!
!! NAME
!!
!!  Hydro_loop1Body
!!
!!
!! SYNOPSIS
!!
!!  Hydro_loop1Body(integer(IN) :: blockCount, 
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
!!  block numbers of blocks on which to advance.
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

!!REORDER(4): scrch_Ptr, scrchFace[XYZ]Ptr, fl[xyz]

Subroutine Hydro_loop1Body(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld,sweepOrder)

  use Eos_interface, ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata,   ONLY : block_metadata_t
  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_putGravityUnsplit

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
                         hy_fullSpecMsFluxHandling

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real, pointer, dimension(:,:,:,:) :: Uout
  real, pointer, dimension(:,:,:,:) :: Uin

  real,dimension(MDIM),intent(IN) :: del
  integer,dimension(LOW:HIGH,MDIM),intent(IN) ::blkLimits,blkLimitsGC 
  type(block_metadata_t), intent(IN) :: blockDesc
  
  integer, dimension(MDIM) :: datasize

  real, allocatable, dimension(:,:,:,:)   :: flx,fly,flz
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)

  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig

  integer :: updateMode ! will be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND


  call Timers_start("loop1 body")
  
  ! NONSENSE...
!!$  call hy_uhd_unsplit(block, Uin,blkLimitsGC,&
!!$                      Uout,blkLimits,&
!!$                      del,dt, dtOld )

!!$     if (hy_fluxCorrect .AND. updateEarly) then
!!$        ! Test whether neighbors are at different refinement levels, and if so,
!!$        ! decide on what we have to do with this block for flux correction.
!!$        call Grid_getBlkNeighLevels(blockID,neighLev)
!!$        myRefine = neighLev(0,0,0)
!!$        blockNeedsFluxCorrect(blockID) = (neighLev(-1,0,0) < myRefine .OR. neighLev(1,0,0) < myRefine)
!!$        blockMustStoreFluxes (blockID) = (neighLev(-1,0,0) .NE. myRefine .OR. neighLev(1,0,0) .NE. myRefine)
!!$#if NDIM > 1
!!$        blockNeedsFluxCorrect(blockID) = blockNeedsFluxCorrect(blockID) &
!!$             .OR. (neighLev(0,-1,0) < myRefine .OR. neighLev(0,1,0) < myRefine)
!!$        blockMustStoreFluxes (blockID) = blockMustStoreFluxes(blockID) &
!!$             .OR. (neighLev(0,-1,0) .NE. myRefine .OR. neighLev(0,1,0) .NE. myRefine)
!!$#endif
!!$#if NDIM > 2
!!$        blockNeedsFluxCorrect(blockID) = blockNeedsFluxCorrect(blockID) &
!!$             .OR. (neighLev(0,0,-1) < myRefine .OR. neighLev(0,0,1) < myRefine)
!!$        blockMustStoreFluxes (blockID) = blockMustStoreFluxes(blockID) &
!!$             .OR. (neighLev(0,0,-1) .NE. myRefine .OR. neighLev(0,0,1) .NE. myRefine)
!!$#endif
!!$     else if (hy_fluxCorrect) then
!!$        blockNeedsFluxCorrect(blockID) = .TRUE.
!!$        blockMustStoreFluxes (blockID) = .TRUE.
!!$     else
!!$        blockNeedsFluxCorrect(blockID) = .FALSE.
!!$        blockMustStoreFluxes (blockID) = .FALSE.
!!$     end if

!!$     call Grid_getDeltas(blockID,del)
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

!!$     if (hy_fullRiemannStateArrays) then
!!$        call hy_memGetBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
!!$        if (NDIM > 1) call hy_memGetBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
!!$        if (NDIM > 2) call hy_memGetBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
!!$     else
!!$     call hy_memAllocScratch(SCRATCH_CTR,HY_VAR1_SCRATCHCTR_VAR,2, 0,0,0, &
!!$          blockList(1:blockCount) )
     allocate(scrch_Ptr(2,dataSize(IAXIS)-1,dataSize(JAXIS)-K2D,dataSize(KAXIS)-K3D))
     allocate(scrchFaceXPtr(HY_NSCRATCH_VARS,dataSize(IAXIS)-1,dataSize(JAXIS)-K2D,dataSize(KAXIS)-K3D))
     allocate(scrchFaceYPtr(HY_NSCRATCH_VARS,dataSize(IAXIS)-1,dataSize(JAXIS)-K2D,dataSize(KAXIS)-K3D))
     allocate(scrchFaceZPtr(HY_NSCRATCH_VARS,dataSize(IAXIS)-1,dataSize(JAXIS)-K2D,dataSize(KAXIS)-K3D))
!!$     endif

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        allocate(  hy_SpcR(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
        allocate(  hy_SpcL(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
        allocate(hy_SpcSig(HY_NSPEC,blkLimits(LOW,IAXIS)-2:blkLimits(HIGH,IAXIS)+2,dataSize(JAXIS),dataSize(KAXIS),NDIM))
     end if
#endif

     allocate(gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#ifdef DEBUG_UHD
     print*,'came upto this point'
#endif
     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        call hy_uhd_putGravityUnsplit(blockDesc,blkLimitsGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ)
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
        call hy_uhd_getRiemannState(blockDesc,Uin,blkLimits,blkLimitsGC,dt,del, &
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
        !    source terms at this point, and will include them in hy_uhd_addGravityUnsplit later
        !    to the primitive Riemann variables (not available for conservative formulation).
        ! 2. With gravity extrapolated from n-1 & n states, gravity source terms have been
        !    included to Riemann states in conservative formulation in hy_uhd_getRiemannState.

     endif !! End of if (hy_updateHydroFluxes) then

     allocate(flx(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NFLUXES))
     allocate(fly(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NFLUXES))
     allocate(flz(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NFLUXES))
     allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))

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
     call hy_uhd_getFaceFlux(blockDesc,blkLimits,blkLimitsGC,datasize,del,&
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
!!$#ifndef FLASH_GRID_UG
!!$     if ((.not. hy_fullRiemannStateArrays) .OR. &
!!$          (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) .OR. &
!!$          .not. blockNeedsFluxCorrect(blockID)) then
!!$#endif
!!$        if (blockNeedsFluxCorrect(blockID)) then
!!$           updateMode = UPDATE_INTERIOR
!!$        else
!!$           updateMode = UPDATE_ALL
!!$        end if
     updateMode = UPDATE_ALL
     call Timers_start("unsplitUpdate")
#ifdef DEBUG_UHD
     print*,'and now update'
#endif
     call hy_uhd_unsplitUpdate(blockDesc,Uin,Uout,updateMode,dt,del,datasize,blkLimits,&
          blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
          scrch_Ptr)
!!$     call io_writeCheckpoint
!!$     stop
!#define DEBUG_UHD
#ifdef DEBUG_UHD
     print*,'done update'
     print*,'_unsplit Aft "call unsplitUpdate(UPD_ALL)": associated(Uin ) is',associated(Uin )
     print*,'_unsplit Aft "call unsplitUpdate(UPD_ALL)": associated(Uout) is',associated(Uout)
#endif
     call Timers_stop("unsplitUpdate")
!!$#ifdef FLASH_UHD_3T
!!$        call Timers_start("unsplitUpdate 3T")
!!$        call hy_uhd_unsplitUpdateMultiTemp&
!!$             (blockID,updateMode,blkLimits,dataSize,dt,del,flx,fly,flz, scrch_Ptr)
!!$        call Timers_stop("unsplitUpdate 3T")
!!$#endif
!!$#ifndef FLASH_GRID_UG
!!$     endif
!!$#endif
!!$     call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)

!!$     if (.not. blockNeedsFluxCorrect(blockID)) then
!!!!!!#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
     call hy_uhd_energyFix(blockDesc,Uout,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
     
#ifdef DEBUG_UHD
     print*,'_unsplit Aft "call energyFix": associated(Uin ) is',associated(Uin )
     print*,'_unsplit Aft "call energyFix": associated(Uout) is',associated(Uout)
#endif
     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
        !! Convert unit
        call hy_uhd_unitConvert(Uout,blkLimitsGC,BWDCONVERT)
     endif
     
     !#ifndef FLASH_EOS_GAMMA
     !! Call to Eos
#ifdef DEBUG_UHD
     print*,'_unsplit bef Eos_wrapped: associated(Uin ) is',associated(Uin )
     print*,'_unsplit bef Eos_wrapped: associated(Uout) is',associated(Uout)
     print*,'_unsplit bef Eos_wrapped: lbound(Uin ):',lbound(Uin )
     print*,'_unsplit bef Eos_wrapped: ubound(Uin ):',ubound(Uin )
     print*,'_unsplit bef Eos_wrapped: lbound(Uout):',lbound(Uout)
     print*,'_unsplit bef Eos_wrapped: ubound(Uout):',ubound(Uout)
#endif
     call Eos_wrapped(hy_eosModeAfter, blkLimits, Uout,CENTER)
     !#endif
!!!!!!#endif /* ifndef GRAVITY */  
     
!!$     if (blockMustStoreFluxes(blockID)) then
        !! if Flux correction is used.
        !! Flux conservation calls on AMR:
        !! Correct fluxes at each block boundary where coarse and fine
        !! blocks are neighboring each other.
        
!!$        if (hy_geometry /= CARTESIAN) then
!!$           ! we are using consv_fluxes and need to divide by face areas
!!$           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
!!$                                (/1,1,1/), faceAreas, datasize)
!!$
!!$           call Grid_putFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
!!$
!!$           if (NDIM > 1) then
!!$              call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
!!$                                   (/1,1,1/), faceAreas, datasize)
!!$              call Grid_putFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
!!$              if (NDIM > 2) then
!!$                 call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
!!$                                      (/1,1,1/), faceAreas, datasize)
!!$                 call Grid_putFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
!!$              endif
!!$           endif
!!$        else ! Cartesian geometry
!!$     call Grid_putFluxData(blockDesc,IAXIS,flx,datasize)
!!$     if (NDIM > 1) then
!!$        call Grid_putFluxData(blockDesc,JAXIS,fly,datasize)
!!$        if (NDIM > 2) then
!!$           call Grid_putFluxData(blockDesc,KAXIS,flz,datasize)
!!$        endif
!!$     endif
!!$     
     
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
     
!!$     if (hy_fullRiemannStateArrays) then
!!$        call hy_memReleaseBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
!!$        if (NDIM > 1) call hy_memReleaseBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
!!$        if (NDIM > 2) call hy_memReleaseBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
!!$     else
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


End Subroutine Hydro_loop1Body
