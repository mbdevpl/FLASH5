#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_ga_isCorrector(ibd, flag)
  use sm_GenAlpha_data, only: sm_GenAlpha_info
  implicit none

  ! IO
  integer, intent(in)  :: ibd 
  integer, intent(out) :: flag

  if( sm_GenAlpha_info(ibd)%pcflag == SM_PCCORRECTOR ) then
     flag = SM_TRUE
  else
     flag = SM_FALSE
  end if

  return

end subroutine sm_ga_isCorrector

