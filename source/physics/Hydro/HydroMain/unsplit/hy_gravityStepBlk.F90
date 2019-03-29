!!****if* source/physics/Hydro/HydroMain/unsplit/hy_gravityStepBlk
!!
!!
!! NAME
!!
!!  hy_gravityStepBlk
!!
!!
!! SYNOPSIS
!!
!!  hy_gravityStepBlk(integer(IN) :: blockCount, 
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld)
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
!!
!!***


Subroutine hy_gravityStepBlk(tileDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld)

  use Eos_interface,    ONLY : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_tile,        ONLY : Grid_tile_t
  use hy_interface,     ONLY : hy_getRiemannState,  &
                               hy_getFaceFlux,      &
                               hy_unsplitUpdate,    &
                               hy_unitConvert,      &
                               hy_energyFix,        &
                               hy_putGravity,&
                               hy_addGravity

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

  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real, pointer, dimension(:,:,:,:) :: Uout
  real, pointer, dimension(:,:,:,:) :: Uin

  real,dimension(MDIM),intent(IN) :: del
  integer,dimension(LOW:HIGH,MDIM),intent(INoUt) ::blkLimits,blkLimitsGC 
  integer :: loxGC,hixGC,loyGC,hiyGC,lozGC,hizGC
  type(Grid_tile_t), intent(IN) :: tileDesc
  
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ

  call Timers_start("loop5 body")


     loxGC = blkLimitsGC(LOW,IAXIS); hixGC =blkLimitsGC(HIGH,IAXIS)
     loyGC = blkLimitsGC(LOW,JAXIS); hiyGC =blkLimitsGC(HIGH,JAXIS)
     lozGC = blkLimitsGC(LOW,KAXIS); hizGC =blkLimitsGC(HIGH,KAXIS)

     allocate(gravX(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravY(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
     allocate(gravZ(loxGC:hixGC, loyGC:hiyGC, lozGC:hizGC))
#ifdef DEBUG_UHD
     print*,'came upto this point in loop5'
#endif
     !! *********************************************************************
     !! Get and add gravity to hydro variables (momenta and energy)         *
     !! *********************************************************************
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        call hy_putGravity(tileDesc,blkLimitsGC,Uin,dt,dtOld,gravX,gravY,gravZ,&
             lastCall=.TRUE.)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref

        call hy_addGravity(tileDesc,blkLimits,blkLimitsGC(LOW,:),blkLimitsGC(HIGH,:),dt,&
             gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:))
     endif


     !! *********************************************************************
     !! Correct energy if necessary
     call hy_energyFix(tileDesc,Uout,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
     
#ifdef DEBUG_UHD
     print*,'_l5 Aft "call energyFix": associated(Uin ) is',associated(Uin )
     print*,'_l5 Aft "call energyFix": associated(Uout) is',associated(Uout)
#endif
     !! *********************************************************************
     !! Correct energy if necessary
     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
        call hy_unitConvert(Uout,blkLimitsGC,BWDCONVERT)
     endif
     
     !! Call to Eos
#ifdef DEBUG_UHD
     print*,'_l5 bef Eos_wrapped: associated(Uin ) is',associated(Uin )
     print*,'_l5 bef Eos_wrapped: associated(Uout) is',associated(Uout)
     print*,'_l5 bef Eos_wrapped: lbound(Uin ):',lbound(Uin )
     print*,'_l5 bef Eos_wrapped: ubound(Uin ):',ubound(Uin )
     print*,'_l5 bef Eos_wrapped: lbound(Uout):',lbound(Uout)
     print*,'_l5 bef Eos_wrapped: ubound(Uout):',ubound(Uout)
#endif
     call Eos_wrapped(hy_eosModeAfter, blkLimits, Uout,CENTER)
     
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
  call Timers_stop("loop5 body")


End Subroutine hy_gravityStepBlk
