!     Stub Function
! File:  
! Author: tim

subroutine sm_assemble_IntForce_rbc(ibd, flag)

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd ! body number
  integer, intent(in)    :: flag
  
  call Driver_abortFlash("sm_Assemble_IntForce_rbc stub function.  Check your Configs")
  
  return
    
end subroutine sm_assemble_IntForce_rbc
