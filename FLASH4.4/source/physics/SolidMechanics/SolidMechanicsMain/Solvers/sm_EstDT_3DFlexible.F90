! sm_pc_EstDT_3DFlexible.F90
!
! Search over each element of a 3DFlexible body and estimate the DT 
! by approximating the largest eigenvalue of each element.

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_EstDT_3DFlexible(ibd, dt)

  use SolidMechanics_data,  only: sm_structure, sm_BodyInfo
  use sm_element_interface, only: sm_3DFlexible_getElement_EvalMax
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,     only: Driver_abortFlash
  implicit none

  ! IO Variables
  integer, intent(in)  :: ibd
  real,    intent(out) :: dt

  ! Internal Variables
  type(sm_structure), pointer :: body
  real    :: Eval_e, Eval_max, TN
  integer :: e

  dt = -1. ! Infamous initialization/

#if NDIM == MDIM

  body => sm_BodyInfo(ibd)

  select case(sm_BodyInfo(ibd)%IntegMethod)
  
  ! Generalized Alpha method
  case(SOLIDINTEG_GENALPHA)
     call RuntimeParameters_get("sminteg3dflexdt",dt)
    
  ! Predictor-Corrector method    
  case(SOLIDINTEG_PREDCORR)
     Eval_max = -9999.

     do e = 1,body%nel
        call sm_3DFlexible_getElement_EvalMax(Eval_e, body, e)
        if( Eval_e > Eval_max ) then
           Eval_max = Eval_e
        end if
     end do

     ! Compute estimate of minimum period
     TN = 2.*PI/sqrt( Eval_max )

     ! Compute conservative DT
     dt = (TN/PI)*0.3
  
  case default
     call Driver_abortFlash('SolidMechanics: Unknown SM dt check scheme.')
     
  end select
       
#endif

  return

end subroutine sm_EstDT_3DFlexible

