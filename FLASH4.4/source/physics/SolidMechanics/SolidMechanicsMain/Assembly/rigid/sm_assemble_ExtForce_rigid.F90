!     Stub Function
! File:  
! Author: tim

subroutine sm_assemble_ExtForce_rigid(ibd, time)

  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "constants.h" 
#include "SolidMechanics.h"   
  !IO Variables
  integer, intent(in)    :: ibd  ! body number
  real, intent(in)       :: time

  ! Local Variables:
  type(sm_structure),  pointer :: body
  integer :: imaster,ix,ex,i,idx,i_dim

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Master node:
  imaster = body % Borigin_node

  if (body%gravity_flag .eq. SM_TRUE) then

     ! ix, ex:
     ix = body % ix; ex = body % ex;

     ! Forces:
     i_dim = 0
     do i = ix,ex
       i_dim = i_dim+1
       idx = body%ID(i,imaster)
       if (idx .le. Body%neq) then
          body%Hs(idx)  = body%Hs(idx) + body%mass * body%gravity_vec(i_dim)
       endif
     end do

  else

     if (ibd .eq. 1) write(*,*) "sm_Assemble_ExtForce_rigid: gravity_flag is .false., no gravity.."  

  endif

!  ! Hack for testing:
!  body%Hs = 0.
!  body%Hs(body%ID(body%ew,imaster)) = .2*sin(2.*PI*1.*time);

  return
    
end subroutine sm_assemble_ExtForce_rigid
