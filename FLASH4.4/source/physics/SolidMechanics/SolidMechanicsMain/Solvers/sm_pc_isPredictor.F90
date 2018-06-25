! Stub

#include "SolidMechanics.h"

subroutine sm_pc_isPredictor(ibd,flag)
  use Driver_interface, only : Driver_abortFlash
  implicit none
  integer, intent(in) :: ibd
  integer, intent(out) :: flag

  call Driver_AbortFlash("PredCorr isPredictor")
  flag = SM_FALSE

end subroutine sm_pc_isPredictor

