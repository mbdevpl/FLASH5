#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_pc_isPredictor(ibd, flag)
  use sm_integdata, only: sm_integ_subiter
  implicit none

  ! IO
  integer, intent(in)  :: ibd 
  integer, intent(out) :: flag

  if( sm_integ_subiter == 0 ) then
     flag = SM_TRUE
  else
     flag = SM_FALSE
  end if

  return

end subroutine sm_pc_isPredictor

