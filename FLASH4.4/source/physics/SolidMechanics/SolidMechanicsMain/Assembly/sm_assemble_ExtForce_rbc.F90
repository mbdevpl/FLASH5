!     Stub Function
! File:  
! Author: tim

subroutine sm_assemble_ExtForce_rbc(ibd, time)

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd! body number
  real, intent(in)       :: time
  
  !call Driver_abortFlash("sm_Assemble_ExtForce_rbc stub function.  Check your Configs")
  
  return
    
end subroutine sm_assemble_ExtForce_rbc
