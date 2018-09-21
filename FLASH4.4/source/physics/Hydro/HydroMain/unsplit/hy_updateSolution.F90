!!****if* source/physics/Hydro/HydroMain/unsplit/hy_updateSolution
!!
!!
!! NAME
!!
!!  hy_updateSolution
!!
!!
!! SYNOPSIS
!!
!!  call hy_updateSolution(block_metadata_t(IN) :: blockDesc,
!!                        integer(IN)          :: blkLimitsGC(LOW:HIGH,MDIM),
!!       real,POINTER(in),dimension(:,:,:,:)   :: Uin,
!!                        integer(IN)          :: blkLimits(LOW:HIGH,MDIM),
!!       real,POINTER(in),dimension(:,:,:,:)   :: Uout,
!!                        real(IN)             :: del(MDM),
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Performs physics update to the solution on a block in a
!!  directionally unsplit fashion.
!!
!!  Fluxes must have been computed before this routine is called.
!!  This implementation expects fluxes to be stored under the
!!  responsibility of the Grid unit, and accessed with pairs of
!!  Grid_getFluxPtr/Grid_releaseFluxPtr calls.
!!
!!  The blockDesc argument tells this routine on which block to
!!  operate.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  blockDesc  - indentifies and describes the current block
!!  Uin        - pointer to one block's worth of cell-centered solution
!!               data; this represents the input data to the current
!!               hydro time step.
!!               Uout must be an associated pointer.
!!               See NOTES below for more info on the relation between
!!               Uin and Uout.
!!  Uout       - pointer to one block's worth of cell-centered solution
!!               data; this represents the output data of the current
!!               hydro time step.
!!               Uout must be an associated pointer.
!!               See NOTES below for more info on the relation between
!!               Uin and Uout.
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a top-layer stub function
!!
!! NOTES
!!
!!  Uin and Uout
!!  ------------
!!
!!  1. Uin and Uout can point to the same storage. Actually this is
!!  the way mostly tested.
!!
!!  2. Schematically (omitting any bothersome details about
!!  directionality, about whether solution variables are in
!!  primitive or conservative form and fluxes are really fluxes or
!!  flux densities, etc.), the function of this routine is to do the
!!  following, in order:
!!
!!  (a) Conservative Update:
!!      **For variables that have corresponding fluxes**
!!
!!        Uout := Uin + fluxes * dt/dx
!!
!!      (not, for example, Uout += fluxes * dt/dx).
!!      **Variables in Uout that do not have corresponding fluxes**,
!!      on the other hand, are not modified in this part (a).
!!      The list of variables that are considered to have
!!      corresponding fluxes is: !!DEV: TBD
!!
!!  (b) Additional Solution Modifications:
!!      These are applied to the variables now in Uout, and can include
!!      * calling hy_unsplitUpdateMultiTemp (for 3T Hydro, not currently
!!        implemented here);
!!      * calling hy_energyFix (to update internal energies, and
!!        possibly other fixups);
!!      * calling hy_unitConvert (for obscure MHD purposes);
!!      * call Eos_wrapped (with mode=hy_eosModeAfter, to bring
!!        solution variables in Uout into a thermodynamically
!!        consistent state).
!!
!!
!!  Other Notes
!!  -----------
!!  The preprocessor symbols MDIM, LOW, HIGH are defined in constants.h .
!!
!!***

!!REORDER(4): scrchFace[XYZ]Ptr

Subroutine hy_updateSolution(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld,sweepOrder)

  use Eos_interface, ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata,   ONLY : block_metadata_t
  use hy_interface, ONLY : hy_getRiemannState,  &
                               hy_getFaceFlux,      &
                               hy_unsplitUpdate,    &
                               hy_unitConvert,      &
                               hy_energyFix,        &
                               hy_putGravity
  use hy_memInterface, ONLY :  hy_memGetBlkPtr,         &
                               hy_memReleaseBlkPtr

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

  use Grid_interface, ONLY : Grid_putFluxData, Grid_getFluxPtr, Grid_releaseFluxPtr
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
  integer :: loxGC,hixGC,loyGC,hiyGC,lozGC,hizGC
  integer :: loFl(MDIM+1)
  type(block_metadata_t), intent(IN) :: blockDesc

  integer, dimension(MDIM) :: datasize

  real, pointer, dimension(:,:,:,:)   :: flx => null() 
  real, pointer, dimension(:,:,:,:)   :: fly => null() 
  real, pointer, dimension(:,:,:,:)   :: flz => null()
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)

  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr => null()
  real, pointer, dimension(:,:,:,:) :: scrchFaceYPtr => null()
  real, pointer, dimension(:,:,:,:) :: scrchFaceZPtr => null()
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr     => null()
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR     => null()
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcL     => null()
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcSig   => null()

  integer :: updateMode ! could be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND



  call Timers_start("update solution body")

  loxGC = blkLimitsGC(LOW,IAXIS); hixGC =blkLimitsGC(HIGH,IAXIS)
  loyGC = blkLimitsGC(LOW,JAXIS); hiyGC =blkLimitsGC(HIGH,JAXIS)
  lozGC = blkLimitsGC(LOW,KAXIS); hizGC =blkLimitsGC(HIGH,KAXIS)
  
#if defined(GPRO_VAR)||defined(VOLX_VAR)||defined(VOLY_VAR)||defined(VOLZ_VAR)||defined(CFL_VAR)
     if (hy_updateHydroFluxes) then
!!$#ifdef GPRO_VAR
!!$        ! A tagging variable for Gaussian Process (GP) method.
!!$        Uin(GPRO_VAR,:,:,:) = 0.
!!$#endif
!!$        !! -----------------------------------------------------------------------!
!!$        !! Save old velocities ---------------------------------------------------!
!!$        !! -----------------------------------------------------------------------!
!!$#ifdef VOLX_VAR
!!$        Uin(VOLX_VAR,:,:,:) = Uin(VELX_VAR,:,:,:)
!!$#endif
!!$#ifdef VOLY_VAR
!!$        Uin(VOLY_VAR,:,:,:) = Uin(VELY_VAR,:,:,:)
!!$#endif
!!$#ifdef VOLZ_VAR
!!$        Uin(VOLZ_VAR,:,:,:) = Uin(VELZ_VAR,:,:,:)
!!$#endif
!!$#ifdef CFL_VAR
!!$        where (1.2*Uin(CFL_VAR,:,:,:) < hy_cfl_original)
!!$           !! Slow recover (of factor of 1.2) to the original CFL once it gets to
!!$           !! reduced to a smaller one in the presence of strong shocks.
!!$           !! This variable CFL takes place in the following three cases using:
!!$           !! (1) use_hybridOrder = .true.,
!!$           !! (2) use_hybridOrder = .true., or
!!$           !! (3) BDRY_VAR is defined and used for stationary objects.
!!$           Uin(CFL_VAR,:,:,:) = 1.2*U(CFL_VAR,:,:,:)
!!$        elsewhere
!!$           Uin(CFL_VAR,:,:,:) = hy_cfl_original
!!$        end where
!!$#endif
!!$     end if
#endif

!!$     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
!!$        call hy_unitConvert(Uin,blkLimitsGC,FWDCONVERT)
!!$     endif

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     call hy_memGetBlkPtr(blockDesc,scrch_Ptr,SCRATCH_CTR)

!!$     allocate(scrchFaceXPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
!!$     allocate(scrchFaceYPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
!!$     allocate(scrchFaceZPtr(HY_NSCRATCH_VARS,loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
!!$
!!$#if (NSPECIES+NMASS_SCALARS) > 0
!!$     if (hy_fullSpecMsFluxHandling) then
!!$        allocate(  hy_SpcR(HY_NSPEC,                                   loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC,NDIM))
!!$        allocate(  hy_SpcL(HY_NSPEC,                                   loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC,NDIM))
!!$        allocate(hy_SpcSig(HY_NSPEC,blkLimits(LOW,IAXIS)-2:blkLimits(HIGH,IAXIS)+2, loyGC:hiyGC, lozGC:hizGC,NDIM))
!!$     end if
!!$#endif
!!$
     allocate(gravX(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravY(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravZ(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
!!$#ifdef DEBUG
!!$     print*,'came upto this point'
!!$#endif
!!$     !! ************************************************************************
!!$     !! Get gravity
!!$     gravX = 0.
!!$     gravY = 0.
!!$     gravZ = 0.
!!$     if (hy_useGravity) then
!!$        call hy_putGravity(blockDesc,blkLimitsGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ)
!!$        gravX = gravX/hy_gref
!!$        gravY = gravY/hy_gref
!!$        gravZ = gravZ/hy_gref
!!$     endif


!!$     if (hy_updateHydroFluxes) then
!!$        !! ************************************************************************
!!$        !! Calculate Riemann (interface) states
!!$        !! Note: gravX(:,:,:) - gravity at n
!!$
!!$#if (NSPECIES+NMASS_SCALARS) > 0
!!$        if (hy_fullSpecMsFluxHandling) then
!!$           hy_SpcL=0.
!!$           hy_SpcR=0.
!!$           hy_SpcSig=0.
!!$        end if
!!$#endif

#ifdef DEBUG_UHD
        print*,'_unsplit bef "call getRiemannState": associated(Uin ) is',associated(Uin )
        print*,'_unsplit bef "call getRiemannState": associated(Uout) is',associated(Uout)
        print*,'_unsplit bef "call getRiemannState": lbound(Uin ):',lbound(Uin )
        print*,'_unsplit bef "call getRiemannState": ubound(Uin ):',ubound(Uin )
#endif
     allocate(  faceAreas(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))

     call Grid_getFluxPtr(blockDesc,flx,fly,flz)
     loFl = lbound(flx)

!!$     if(hy_fluxCorrectPerLevel) then
!!$        updateMode=UPDATE_ALL
!!$        call Grid_getFluxData(blockDesc,flx,fly,flz,datasize)
!!$     else
!!$        
!!$        updateMode = UPDATE_INTERIOR
!!$     end if

     updateMode = UPDATE_ALL

     call Timers_start("unsplitUpdate")
#ifdef DEBUG_UHD
     print*,'and now update'
#endif
     
     call hy_unsplitUpdate(blockDesc,Uin,Uout,updateMode,dt,del,datasize,blkLimits,&
          blkLimitsGC,loFl,flx,fly,flz,gravX,gravY,gravZ,&
          scrch_Ptr)
     
     
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
     call hy_memReleaseBlkPtr(blockDesc,scrch_Ptr,SCRATCH_CTR)

!!$     if (.not. blockNeedsFluxCorrect(blockID)) then
#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
     call hy_energyFix(blockDesc,Uout,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
     
#ifdef DEBUG_UHD
     print*,'_unsplit Aft "call energyFix": associated(Uin ) is',associated(Uin )
     print*,'_unsplit Aft "call energyFix": associated(Uout) is',associated(Uout)
#endif
     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
        !! Convert unit
        call hy_unitConvert(Uout,blkLimitsGC,BWDCONVERT)
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
#endif /* ifndef GRAVITY */
     
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

!!$     if (hy_fluxCorrect) then
!!$        call Grid_putFluxData(blockDesc,flx,fly,flz,datasize)
!!$     end if
     call Grid_releaseFluxPtr(blockDesc,flx,fly,flz)
     
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
     
  call Timers_stop("update solution body")


End Subroutine hy_updateSolution
