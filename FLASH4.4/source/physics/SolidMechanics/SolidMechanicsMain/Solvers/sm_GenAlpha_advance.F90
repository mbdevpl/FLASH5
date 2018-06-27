! Stub

subroutine sm_GenAlpha_advance(ibd,restart)
  use Driver_interface, only : Driver_abortFlash
  implicit none
  integer, intent(in) :: ibd
  logical, intent(in) :: restart

  call Driver_AbortFlash("GenAlpha_init not configured")

end subroutine sm_GenAlpha_advance

