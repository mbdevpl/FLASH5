! Stub

subroutine sm_ga_isCorrector(ibd,flag)
  use Driver_interface, only : Driver_abortFlash
  implicit none
  integer, intent(in) :: ibd
  integer, intent(out) :: flag

  call Driver_AbortFlash("GenAlpha isCorrector")
  flag = -1

end subroutine sm_ga_isCorrector

