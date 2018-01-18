!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/hy_uhd_unsplit
!!
!! NAME
!!
!!  hy_uhd_unsplit
!!
!! SYNOPSIS
!!
!!  call hy_uhd_unsplit(integer (IN) :: blockCount,
!!                  integer (IN) :: blockList(blockCount),
!!                  real    (IN) :: dt,
!!                  real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs Hydro update in a directionally unsplit fashion over a set
!!  of blocks.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - Eos_wrapped call is applied to the guard cells where necessary;
!!   - computes fluxes using a call to hy_uhd_getFaceFlux
!!   - if we're not doing flux correction (as controlled by the flux_correct
!!     runtime parameter), then we update all the cell values from the fluxes 
!!     (with a call to hy_uhd_unsplitUpdate), otherwise, we update just cells 
!!     not on the boundaries, and save fluxes for cells on the boundary;
!!   - and finally, we apply an eos to the block.
!! 
!!  After the main block loop, if doing flux correction, we have
!!  the Grid correct boundary fluxes for all blocks where approriate,
!!  and do another loop over blocks, updating the cell values for
!!  cells on the block boundaries using the corrected fluxes, and
!!  apply an eos on the block. 
!!
!!
!! REFERENCES
!!
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Lee, D., "A solution accurate, efficient and stable unsplit staggered mesh scheme 
!!              for three dimensional magnetohydrodynamics", 243 (2013), 269-292, JCP
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (needed by some implementations for temporal extrapolations
!!                              of gravity; not used in this implementation)
!!
!!***

!!REORDER(4): U, scrch_Ptr, scrchFace[XYZ]Ptr, fl[xyz]

#ifdef DEBUG_ALL
#define DEBUG_UHD
#endif
#define DEBUG_GRID_GCMASK

#include "Flash.h"
Subroutine hy_uhd_unsplit (blockDesc,Uin,blkLimitsGC,&
                      Uout,blkLimits,&
                      del,dt, dtOld )

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


  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime

  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_prepareNewGravityAccel,&
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_shockDetect,      &
                               hy_uhd_unsplitUpdateMultiTemp,&
                               hy_uhd_multiTempAfter
  use hy_memInterface, ONLY :  hy_memAllocScratch,      &
                               hy_memDeallocScratch,    &
                               hy_memGetBlkPtr,         &
                               hy_memReleaseBlkPtr

  use Grid_interface, ONLY : Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_getBlkData,        &
                             Grid_getBlkNeighLevels

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata,   ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  type(block_metadata_t), intent(IN) :: blockDesc
  real,    INTENT(IN) :: dt, dtOld
  integer,dimension(LOW:HIGH,MDIM),INTENT(IN) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:),pointer :: Uin
  real,dimension(:,:,:,:), pointer :: Uout
  real,dimension(MDIM),INTENT(IN) :: del
 !! -----------------------------------------------------

  integer, dimension(MDIM) :: datasize
  integer, parameter :: level=0
  integer :: i,blockID
  integer :: ix,iy,iz
  logical :: gcMask(hy_gcMaskSize)

  real, allocatable, dimension(:,:,:,:)   :: flx,fly,flz
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)


  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  real, pointer, dimension(:,:,:,:) :: U

  integer :: updateMode ! will be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND
  logical,dimension(MAXBLOCKS) :: blockNeedsFluxCorrect, blockMustStoreFluxes
  logical, parameter :: updateEarly=.FALSE. !update in the first loop if possible?
  integer :: neighLev(-1:1, -K2D:K2D , -K3D:K3D)
  integer :: myRefine
  !! End of data declaration ***********************************************

#ifdef DEBUG_UHD
  print*,'_unsplit entry: associated(Uin ) is',associated(Uin )
  print*,'_unsplit entry: associated(Uout) is',associated(Uout)
#endif


  
  
  
  
! ********** STUFF MOVED TO Hydro_prepareBuffers **********

#ifdef DEBUG_UHD
  print*,'_unsplit bef "LOOP 1": associated(Uin ) is',associated(Uin )
  print*,'_unsplit bef "LOOP 1": associated(Uout) is',associated(Uout)
#endif


!!$  do i=1,blockCount             !LOOP 1
!!$
!!$     blockID = blockList(i)


! ********** STUFF MOVED TO Hydro_loop1Body **********

     
     
     
!!$     !! ***************************************************************************
!!$     !! Third part of advancement                                                 *
!!$     !! ***************************************************************************
!!$     !! Do this part only if refining and flux correcting
!!$     if (hy_fluxCorrect) then
!!$        
!!$        !! ************************************************************************
!!$        !! Conservation of Fluxes at each block boundary
!!$        call Timers_start("conserveFluxes")
!!$        call Grid_conserveFluxes(ALLDIR,level)
!!$        call Timers_stop("conserveFluxes")
!!$        
!!$
!!$        datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
!!$        allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(    gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(    gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(    gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$        allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$
!!$        !! *********************************************************************
!!$        !! Get gravity
!!$        gravX = 0.
!!$        gravY = 0.
!!$        gravZ = 0.
!!$        if (hy_useGravity) then
!!$           call hy_uhd_putGravityUnsplit(blockDesc,Uout,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
!!$           gravX = gravX/hy_gref
!!$           gravY = gravY/hy_gref
!!$           gravZ = gravZ/hy_gref
!!$        endif
!!$
!!$
!!$        !! Get modified conserved flux values at block interfaces:
!!$        !! This is important especially at the block interface at 
!!$        !! different refinement levels of neighboring blocks
!!$        flx = 0.
!!$        fly = 0.
!!$        flz = 0.
!!$
!!$        if (hy_fullRiemannStateArrays) then
!!$           call hy_memGetBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
!!$           if (NDIM > 1) call hy_memGetBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
!!$           if (NDIM > 2) call hy_memGetBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
!!$        endif
!!$
!!$        if (hy_fullRiemannStateArrays) then
!!$           flx(HY_DENS_FLUX:HY_END_FLUX,&
!!$                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
!!$                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
!!$              =scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,&
!!$                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
!!$                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
!!$#if NDIM > 1
!!$           fly(HY_DENS_FLUX:HY_END_FLUX,&
!!$                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
!!$                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
!!$              =scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,&
!!$                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
!!$                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
!!$#if NDIM == 3
!!$           flz(HY_DENS_FLUX:HY_END_FLUX,&
!!$                blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)&
!!$              =scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,&
!!$                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
!!$#endif
!!$#endif
!!$        endif
!!$
!!$        if (hy_geometry /= CARTESIAN) then
!!$           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
!!$                                (/1,1,1/), faceAreas, datasize)
!!$
!!$           call Grid_getFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
!!$
!!$        else
!!$           call Grid_getFluxData(blockDesc,IAXIS,flx,datasize)
!!$        endif
!!$
!!$#if NDIM > 1
!!$        if (hy_geometry /= CARTESIAN) then
!!$           call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
!!$                                (/1,1,1/), faceAreas, datasize)
!!$           call Grid_getFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
!!$        else
!!$           call Grid_getFluxData(blockDesc,JAXIS,fly,datasize)
!!$        endif
!!$#if NDIM > 2
!!$        if (hy_geometry /= CARTESIAN) then
!!$           call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
!!$                                (/1,1,1/), faceAreas, datasize)
!!$           call Grid_getFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
!!$        else
!!$           call Grid_getFluxData(blockID,KAXIS,flz,datasize)
!!$        endif
!!$#endif
!!$#endif
!!$
!!$        if (hy_fullRiemannStateArrays) then
!!$           call hy_memReleaseBlkPtr(blockID,scrchFaceXPtr,SCRATCH_FACEX)
!!$           if (NDIM > 1) call hy_memReleaseBlkPtr(blockID,scrchFaceYPtr,SCRATCH_FACEY)
!!$           if (NDIM > 2) call hy_memReleaseBlkPtr(blockID,scrchFaceZPtr,SCRATCH_FACEZ)
!!$        endif
!!$
!!$        call hy_memGetBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
!!$        !! *********************************************************************
!!$        !! Unsplit update for conservative variables from n to n+1 time step
!!$        if (.not. hy_fullRiemannStateArrays .OR. (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0)) then
!!$           updateMode=UPDATE_BOUND
!!$        else
!!$           updateMode=UPDATE_ALL
!!$        endif
!!$        call Timers_start("unsplitUpdate")
!!$        call hy_uhd_unsplitUpdate(blockDesc,updateMode,dt,del,datasize,blkLimits,&
!!$                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ,&
!!$                                  scrch_Ptr)
!!$        call Timers_stop("unsplitUpdate")
!!$#ifdef FLASH_UHD_3T
!!$        call Timers_start("unsplitUpdate 3T")
!!$        call hy_uhd_unsplitUpdateMultiTemp&
!!$             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz, scrch_Ptr)
!!$        call Timers_stop("unsplitUpdate 3T")
!!$#endif
!!$
!!$#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
!!$        !! *********************************************************************
!!$        !! Correct energy if necessary
!!$        call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
!!$
!!$        !! *********************************************************************
!!$        !! Convert units if necessary
!!$        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
!!$           !! Convert units
!!$           call hy_uhd_unitConvert(Uout,blkLimitsGC,BWDCONVERT)
!!$        endif
!!$
!!$
!!$        !! *********************************************************************
!!$        !! Call to Eos
!!$        call Eos_wrapped(hy_eosModeAfter, blkLimits, Uout,CENTER)
!!$#endif /* if gravity is included, we calculate gravity for n+1 in the below */
!!$
!!$        call hy_memReleaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
!!$
!!$        deallocate(flx)
!!$        deallocate(fly)
!!$        deallocate(flz)
!!$        deallocate(gravX)
!!$        deallocate(gravY)
!!$        deallocate(gravZ)
!!$        deallocate(faceAreas)
!!$        
!!$     end if !! End of the flux conservation routine
!!$
!!$     call Timers_start("tail")
!!$     
!!$     if (hy_fullRiemannStateArrays) then
!!$        call hy_memDeallocScratch(SCRATCH_FACEX)
!!$        if (NDIM > 1) call hy_memDeallocScratch(SCRATCH_FACEY)
!!$        if (NDIM > 2) call hy_memDeallocScratch(SCRATCH_FACEZ)
!!$     endif
!!$
!!$
!!$#ifdef GRAVITY /* Perform this only when gravity is used */
!!$  !! ***************************************************************************
!!$  !! Fourth part of advancement to compute gravity at n+1 state                *
!!$  !! ***************************************************************************
!!$
!!$#ifdef GPOT_VAR
!!$  if (hy_useGravity) then
!!$     ! The following call invokes Gravity_potentialListOfBlocks and related stuff,
!!$     ! to prepare for retrieving updated accelerations below.
!!$     call hy_uhd_prepareNewGravityAccel(blockCount,blockList,gcMaskLogged)
!!$  endif
!!$#endif
!!$
!!$
!!$  !! Proceed to couple the updated gravitational accelerations 
!!$  !! to energy and momenta on each block (see hy_uhd_addGravityUnsplit)
!!$  do i = 1,blockCount           !LOOP 5
!!$
!!$     blockID = blockList(i)
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!$     call Grid_getDeltas(blockID,del)
!!$
!!$     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
!!$#ifndef FIXEDBLOCKSIZE
!!$     allocate(gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$     allocate(gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$     allocate(gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
!!$#endif
!!$
!!$     !! *********************************************************************
!!$     !! Get and add gravity to hydro variables (momenta and energy)         *
!!$     !! *********************************************************************
!!$     gravX = 0.
!!$     gravY = 0.
!!$     gravZ = 0.
!!$     if (hy_useGravity) then
!!$
!!$        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ,&
!!$             lastCall=.TRUE.)
!!$        gravX = gravX/hy_gref
!!$        gravY = gravY/hy_gref
!!$        gravZ = gravZ/hy_gref
!!$
!!$        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,dt,&
!!$             gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:))
!!$     endif
!!$
!!$     !! *********************************************************************
!!$     !! Correct energy if necessary
!!$     call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
!!$
!!$     !! *********************************************************************
!!$     !! Convert units if necessary
!!$     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
!!$        !! Convert units
!!$        call hy_uhd_unitConvert(blockID,BWDCONVERT)
!!$     endif
!!$
!!$     !! *********************************************************************
!!$     !! Call to Eos
!!$     call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
!!$
!!$#ifndef FIXEDBLOCKSIZE
!!$     deallocate(gravX)
!!$     deallocate(gravY)
!!$     deallocate(gravZ)
!!$#endif
!!$
!!$  end do !! End of the loop over blocks
!!$#endif /* End of n+1 gravity coupling */
!!$
!!$
!!$
!!$
!!$  ! Do 3T update...
!!$#ifdef FLASH_UHD_3T
!!$  call hy_uhd_multiTempAfter(blockCount, blockList, dt)
!!$#endif
!!$

! ********** STUFF MOVED TO Hydro_freeBuffers  **********


! ********** STUFF MOVED TO Hydro_advanceAll  **********



  call Timers_stop("tail")

End Subroutine hy_uhd_unsplit
