!     Stub Function
! File:  
! Author: tim

subroutine sm_assemble_stiff_rbc(ibd, flag)

  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd, flag ! body number
  
  call Driver_abortFlash("sm_Assemble_stiff_rbc stub function.  Check your Configs")
  
  return
    
end subroutine sm_assemble_stiff_rbc
