!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_putGravityUnsplit
!!
!! NAME
!!
!!  hy_uhd_PutGravityUnsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_putGravityUnsplit( integer(IN) :: Uin,
!!                            integer(IN) :: blkLimitsGC(2,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real   (IN) :: dt,
!!                            real   (IN) :: dtOld,
!!                            real(OUT)   :: gravX(:,:,:),
!!                            real(OUT)   :: gravY(:,:,:),
!!                            real(OUT)   :: gravZ(:,:,:))
!!
!! ARGUMENTS
!!
!!  blockID     - a current block ID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dataSize    - dimensions for gravX, gravY and gravZ arrays
!!  dt          - timestep
!!  dtOld       - old timestep (needed for temporal extrapolations of gravity)
!!  gravX       - gravity components in x-direcition at time steps n
!!  gravY       - gravity components in y-direcition at time steps n
!!  gravZ       - gravity components in z-direcition at time steps n
!!
!! DESCRIPTION
!!
!!  This routine puts gravity components to arrays gravX, gravY and gravZ:
!!  gravX(:,:,:) includes gravity components at time step n
!!
!!*** 

Subroutine hy_uhd_putGravityUnsplit(blockDesc,blGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ, potentialIndex, lastCall)

  use Gravity_interface, ONLY : Gravity_accelOneRow

  use Hydro_data, ONLY: hy_gpotVar, hy_extraAccelVars
  use block_metadata,   ONLY : block_metadata_t
  
  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  type(block_metadata_t), intent(IN)   :: blockDesc
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blGC
  real,dimension(:,:,:,:),pointer :: Uin
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt, dtOld

  real, dimension(blGC(LOW,IAXIS):blGC(HIGH,IAXIS), blGC(LOW,JAXIS):blGC(HIGH,JAXIS), blGC(LOW,KAXIS):blGC(HIGH,KAXIS)), &
       intent(OUT) :: gravX,gravY,gravZ

  integer, intent(IN), OPTIONAL :: potentialIndex
  logical, intent(IN), OPTIONAL :: lastCall
  !! -----------------------------------------------------

  integer, dimension(2) :: gravPos
  integer :: ix, iy, iz
  integer :: potVar
  logical :: tailCall

  if (present(lastCall)) then
     tailCall = lastCall
  else
     tailCall = .FALSE.
  end if

  do iz=blGC(LOW,KAXIS),blGC(HIGH,KAXIS)
     do iy=blGC(LOW,JAXIS),blGC(HIGH,JAXIS)
        gravPos(1)=iy
        gravPos(2)=iz
        if (tailCall .AND. hy_gpotVar > 0) then
           if (present(potentialIndex)) then
              potVar = potentialIndex
           else
              potVar = hy_gpotVar
           end if
           if (hy_extraAccelVars(1)>0) then
              call Gravity_accelOneRow(gravPos,DIR_X,blockDesc,dataSize(IAXIS),gravX(:,iy,iz),Uin,potVar,&
                                       extraAccelVars=hy_extraAccelVars)
           else
              call Gravity_accelOneRow(gravPos,DIR_X,blockDesc,dataSize(IAXIS),gravX(:,iy,iz),Uin,potVar)
           end if
        else if (present(potentialIndex)) then
           call Gravity_accelOneRow(gravPos,DIR_X,blockDesc,dataSize(IAXIS),gravX(:,iy,iz),Uin,potentialIndex)
        else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
        ! Gravity implementation defines FLASH_GRAVITY_TIMEDEP -> time-dependent gravity field
        ! gravity at time step n
           call Gravity_accelOneRow(gravPos,DIR_X,blockDesc,dataSize(IAXIS),gravX(:,iy,iz),Uin,GPOT_VAR)
#else
        ! FLASH_GRAVITY_TIMEDEP not defined -> assume time-independent gravity field.
        ! Also if GPOT_VAR is defined -> use current accel without time
        ! interpolation, i.e., handle like time-independent gravity field - KW
           call Gravity_accelOneRow(gravPos,DIR_X,blockDesc,dataSize(IAXIS),gravX(:,iy,iz),Uin)
#endif
        endif
     enddo
  enddo


  if (NDIM >= 2) then
     do iz=blGC(LOW,KAXIS),blGC(HIGH,KAXIS)
        do ix=blGC(LOW,IAXIS),blGC(HIGH,IAXIS)
           gravPos(1)=ix
           gravPos(2)=iz
           if (tailCall .AND. hy_gpotVar > 0) then
              if (present(potentialIndex)) then
                 potVar = potentialIndex
              else
                 potVar = hy_gpotVar
              end if
              if (hy_extraAccelVars(2)>0) then
                 call Gravity_accelOneRow(gravPos,DIR_Y,blockDesc,dataSize(JAXIS),gravY(ix,:,iz),Uin,potVar,&
                                          extraAccelVars=hy_extraAccelVars)
              else
                 call Gravity_accelOneRow(gravPos,DIR_Y,blockDesc,dataSize(JAXIS),gravY(ix,:,iz),Uin,potVar)
              end if
           else if (present(potentialIndex)) then
              call Gravity_accelOneRow(gravPos,DIR_Y,blockDesc,dataSize(JAXIS),gravY(ix,:,iz),Uin,potentialIndex)
           else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
           ! gravity at time step n
              call Gravity_accelOneRow(gravPos,DIR_Y,blockDesc,dataSize(JAXIS),gravY(ix,:,iz),Uin,GPOT_VAR)
#else
              call Gravity_accelOneRow(gravPos,DIR_Y,blockDesc,dataSize(JAXIS),gravY(ix,:,iz),Uin)
#endif
           end if
        enddo
     enddo


     if (NDIM == 3) then
        do iy=blGC(LOW,JAXIS),blGC(HIGH,JAXIS)
           do ix=blGC(LOW,IAXIS),blGC(HIGH,IAXIS)
              gravPos(1)=ix
              gravPos(2)=iy
              if (tailCall .AND. hy_gpotVar > 0) then
                 if (present(potentialIndex)) then
                    potVar = potentialIndex
                 else
                    potVar = hy_gpotVar
                 end if
                 if (hy_extraAccelVars(3)>0) then
                    call Gravity_accelOneRow(gravPos,DIR_Z,blockDesc,dataSize(KAXIS),gravZ(ix,iy,:),Uin,potVar,&
                                             extraAccelVars=hy_extraAccelVars)
                 else
                    call Gravity_accelOneRow(gravPos,DIR_Z,blockDesc,dataSize(KAXIS),gravZ(ix,iy,:),Uin,potVar)
                 end if
              else if (present(potentialIndex)) then
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockDesc,dataSize(KAXIS),gravZ(ix,iy,:),Uin,potentialIndex)
              else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
              ! gravity at time step n
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockDesc,dataSize(KAXIS),gravZ(ix,iy,:),Uin,GPOT_VAR)
#else
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockDesc,dataSize(KAXIS),gravZ(ix,iy,:),Uin)
#endif
              end if
           enddo
        enddo
     endif
  endif
End Subroutine hy_uhd_putGravityUnsplit
