!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_integ_adv1dt
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!!
!!  
!!
!!
!! DESCRIPTION
!! Advances one subiteration, or timestep depending on the time integration desired.
!! For Predictor Corrector-schemes does one predictor or one of the correction 
!! subiterations.
!!
!!***
#include "SolidMechanics.h"

subroutine sm_integ_adv1dt(restart_local)

  use SolidMechanics_data, only : sm_MeshMe
  use sm_integinterface, only : sm_integ_advance, sm_integ_checkconverg
  use sm_integdata, only : sm_convflag_all, sm_errmax_all, sm_integ_subiter 
  implicit none
  logical, INTENT(IN) :: restart_local
  
  integer :: convflag_all


  ! Reset Global integration variables for new timestep
  sm_convflag_all = SM_NOTCONVERGED
  convflag_all    = sm_convflag_all
  sm_errmax_all   = 0.
  sm_integ_subiter= 0

  ! Iteration to convergence:
  do while (convflag_all .ne. SM_CONVERGED)

     ! Advance substep:
     call sm_integ_advance(restart_local)


     ! Check convergence:
     call sm_integ_checkconverg(convflag_all)

  end do

  return

end subroutine sm_integ_adv1dt
