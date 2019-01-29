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
!!  call hy_computeFluxes(block_metadata_t(IN) :: blockDesc,
!!                        integer(IN)          :: blkLimitsGC(LOW:HIGH,MDIM),
!!       real,POINTER(in),dimension(:,:,:,:)   :: Uin,
!!                        integer(IN)          :: blkLimits(LOW:HIGH,MDIM),
!!       real,POINTER(in),dimension(:,:,:,:)   :: Uout,
!!                        real(IN)             :: del(MDM),
!!                        real(IN)             :: timeEndAdv,
!!                        real(IN)             :: dt,
!!                        real(IN)             :: dtOld,
!!                        integer(IN)          :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Computes hydrodynamical fluxes for all spatial directions for
!!  a block, and optionally applies those fluxes to update the
!!  cell-centered solution in UNK.
!!
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
!!  blockDesc  - describes the current block
!!  Uin        - pointer to one block's worth of cell-centered solution
!!               data; this represents the input data to the current
!!               hydro time step.
!!  Uout       - pointer to one block's worth of cell-centered solution
!!               data; this represents the output data of the current
!!               hydro time step.
!!               Uout may be a disassociated pointer. In that case,
!!               fluxes are computed but not applied to update the
!!               cell-centered solution.
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a top-layer stub function
!!
!! NOTES
!!
!!  Contrary to what the name might suggests, this code includes the action
!!  of the also separately available routine hy_updateSolution when
!!  the pointer Uout is non-null.
!!
!!  The preprocessor symbols MDIM, LOW, HIGH are defined in constants.h .
!!
!! SEE ALSO
!!
!!  hy_updateSolution
!!
!!***


!!$#define DEBUG
!!$#define DEBUG_UHD

!!REORDER(4): scrch_Ptr, scrchFace[XYZ]Ptr, fl[xyz]

Subroutine hy_computeFluxes(tileDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld,sweepOrder)
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use flash_tile,   ONLY : flash_tile_t 
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

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real, pointer, dimension(:,:,:,:) :: Uout
  real, pointer, dimension(:,:,:,:) :: Uin

  real,dimension(MDIM),intent(IN) :: del
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
  integer :: loxGC,hixGC,loyGC,hiyGC,lozGC,hizGC
  integer :: loFl(MDIM+1)
  type(flash_tile_t), intent(IN) :: tileDesc

  integer, dimension(MDIM) :: datasize

  real, pointer, dimension(:,:,:,:)   :: flx
  real, pointer, dimension(:,:,:,:)   :: fly
  real, pointer, dimension(:,:,:,:)   :: flz
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, target, dimension(0,0,0,0)   :: empty4

  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr
  real, pointer, dimension(:,:,:,:) :: scrchFaceYPtr
  real, pointer, dimension(:,:,:,:) :: scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr    
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR  
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcL  
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcSig

  integer :: updateMode ! will be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND

  integer :: halo(LOW:HIGH, 1:MDIM)

  scrchFaceXPtr => null()
  scrchFaceYPtr => null()
  scrchFaceZPtr => null()
  scrch_Ptr     => null()
  flx           => null() 
  fly           => null()
  flz           => null()
  hy_SpcR       => null() 
  hy_SpcL       => null()
  hy_SpcSig     => null()

  call Timers_start("loop1 body")

  if (associated(Uout)) then
     updateMode = UPDATE_ALL    ! Flag that we should do the update immediately here, too!
  else
     updateMode = UPDATE_NONE
  end if

  ! Temporarily assume that all internal arrays created here and in routines
  ! called here will need data for the tile interior as well as for NGUARD
  ! layers of GCs around the interior
  halo(:, :) = tileDesc%limits(:, :)
  halo(LOW,  1:NDIM) = halo(LOW,  1:NDIM) - NGUARD
  halo(HIGH, 1:NDIM) = halo(HIGH, 1:NDIM) + NGUARD
  loxGC = halo(LOW,IAXIS); hixGC = halo(HIGH,IAXIS)
  loyGC = halo(LOW,JAXIS); hiyGC = halo(HIGH,JAXIS)
  lozGC = halo(LOW,KAXIS); hizGC = halo(HIGH,KAXIS)

#if defined(GPRO_VAR)||defined(VOLX_VAR)||defined(VOLY_VAR)||defined(VOLZ_VAR)||defined(CFL_VAR)
     if (hy_updateHydroFluxes) then
        ! Set with explicit loop nests over a tile-based region
        call Driver_abortFlash("Update this to work with tiles")
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
        call Driver_abortFlash("Confirm that grownLimits is correct")
        call hy_unitConvert(Uin, tileDesc%grownLimits, FWDCONVERT)
     endif

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     if (updateMode == UPDATE_ALL) then
        allocate(scrch_Ptr    (2,            loxGC:hixGC-1, loyGC:hiyGC-K2D, lozGC:hizGC-K3D))
     else
        call hy_memGetBlkPtr(tileDesc,scrch_Ptr,SCRATCH_CTR)
     end if
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
        call Driver_abortFlash("Not using Gravity yet!")
        call hy_putGravity(tileDesc,blkLimitsGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ)
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
        call hy_getRiemannState(tileDesc, Uin, tileDesc%limits, &
                                halo(LOW,:), halo(HIGH,:), &
                                dt, del, &
                                gravX, gravY, gravZ,&
                                scrchFaceXPtr, scrchFaceYPtr, scrchFaceZPtr,&
                                hy_SpcR, hy_SpcL, hy_SpcSig)
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
     call tileDesc%getDataPtr(flx, FLUXX)
     call tileDesc%getDataPtr(fly, FLUXY)
     call tileDesc%getDataPtr(flz, FLUXZ)
     loFl = lbound(flx)
     if (.NOT. associated(fly)) fly => empty4
     if (.NOT. associated(flz)) flz => empty4

!!$     call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR) 

     !! ************************************************************************
     !! Calculate high order Godunov fluxes by making Riemann solver calls
     call Timers_start("getFaceFlux")
#ifdef DEBUG_UHD
     print*,'getting face flux'
#endif
     call hy_getFaceFlux(tileDesc, tileDesc%limits, halo, &
                         datasize, del, &
                         loFl, flx, fly, flz, &
                         scrchFaceXPtr, scrchFaceYPtr, scrchFaceZPtr, &
                         scrch_Ptr, hy_SpcR, hy_SpcL)
#ifdef DEBUG_UHD
     print*,'got face flux'
     print*,'_unsplit Aft "call getFaceFlux": associated(Uin ) is',associated(Uin )
     print*,'_unsplit Aft "call getFaceFlux": associated(Uout) is',associated(Uout)
#endif
     call Timers_stop("getFaceFlux")
     if (updateMode == UPDATE_ALL) then
        call Timers_start("unsplitUpdate")
#ifdef DEBUG_UHD
        print*,'and now update'
#endif
        call hy_unsplitUpdate(tileDesc, Uin, Uout, updateMode, &
                              dt, del, datasize, &
                              tileDesc%limits, halo, &
                              loFl, flx, fly, flz, &
                              gravX, gravY, gravZ, &
                              scrch_Ptr)
        call Timers_stop("unsplitUpdate")
#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
        call hy_energyFix(tileDesc, Uout, tileDesc%limits, dt, dtOld, &
                          del, hy_unsplitEosMode)

#ifdef DEBUG_UHD
        print*,'_unsplit Aft "call energyFix": associated(Uin ) is',associated(Uin )
        print*,'_unsplit Aft "call energyFix": associated(Uout) is',associated(Uout)
#endif
        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call Driver_abortFlash("Confirm that grownLimits is correct")
           call hy_unitConvert(Uout, tileDesc%grownLimits, BWDCONVERT)
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
        call Eos_wrapped(hy_eosModeAfter, tileDesc%limits, Uout,CENTER)
     !#endif
#endif /* ifndef GRAVITY */

        deallocate(scrch_Ptr)
     else
        call hy_memReleaseBlkPtr(tileDesc,scrch_Ptr,SCRATCH_CTR)
     end if
     !! ************************************************************************
     !! Unsplit update for conservative variables from n to n+1 time step

     
     call tileDesc%releaseDataPtr(flx, FLUXX)
     call tileDesc%releaseDataPtr(fly, FLUXY)
     call tileDesc%releaseDataPtr(flz, FLUXZ)
     
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     
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
