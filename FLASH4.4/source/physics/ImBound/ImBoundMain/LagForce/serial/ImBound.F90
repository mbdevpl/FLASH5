!!****if* source/physics/ImBound/ImBoundMain/LagForce/serial/ImBound
!!
!!
!! NAME
!!
!!  ImBound
!!
!!
!! SYNOPSIS
!!
!!  ImBound(blockCount,blockList,dt)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ImBound(blockCount, blockList, dt, forcflag)

  use Driver_interface, only : Driver_abortFlash

  use ib_interface, only : ib_forcing,ib_CalcForce

  implicit none
#include "Flash.h"
#include "ImBound.h"
  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  integer, INTENT(IN) :: forcflag
  !! -----------------------------------------------------

  select case (forcflag)
  case(FORCE_FLOW) 

  call ib_forcing(blockCount, blockList, dt)

  case(COMPUTE_FORCES)

  call ib_CalcForce(blockCount, blockList, dt)

  case default

  call Driver_abortFlash("ImBound : forcflag doen not correspond to any available option.") 
  
  end select

  return
end subroutine ImBound

