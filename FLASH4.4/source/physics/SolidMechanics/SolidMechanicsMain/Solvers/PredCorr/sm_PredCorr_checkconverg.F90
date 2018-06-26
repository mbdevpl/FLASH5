!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/PredCorr/sm_PredCorr_checkconverg
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
!! Checks for convergence on Predictor-Corrector schemes.
!!
!!***
#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_PredCorr_checkconverg(testflg,ibd,convflag_pc,errmax_pc)

  use Driver_interface, only : Driver_abortFlash
  use sm_PredCorr_data, only : sm_PredCorr_Info

  implicit none
  integer, intent(IN)  :: testflg,ibd
  integer, intent(OUT) :: convflag_pc
  real, intent(OUT)    :: errmax_pc


  convflag_pc = -99999
  errmax_pc   =  1.

  select case(testflg)

  case(SM_TESTCONVERGE) 

    convflag_pc = sm_PredCorr_info(ibd)%pcconvflag
    errmax_pc   = sm_PredCorr_info(ibd)%pcerr

  case(SM_SETNOTCONVERGED) ! Case one body sub-iteration has not converged.

    sm_PredCorr_info(ibd)%pcconvflag = SM_PCNOTCONVERGED
    sm_PredCorr_info(ibd)%pcflag     = SM_PCCORRECTOR
 
    ! IF at start up we might need to lower the method in pcmethod !

  case default

     call Driver_abortFlash('sm_PredCorr_checkconverg: Unknown test flag.')
  
  end select

  return

end subroutine sm_PredCorr_checkconverg
