! Stub

subroutine sm_ga_isPredictor(ibd,flag)
  use Driver_interface, only : Driver_abortFlash
  implicit none
  integer, intent(in) :: ibd
  integer, intent(out) :: flag

  call Driver_AbortFlash("GenAlpha isPredictor")
  flag = -1

end subroutine sm_ga_isPredictor

