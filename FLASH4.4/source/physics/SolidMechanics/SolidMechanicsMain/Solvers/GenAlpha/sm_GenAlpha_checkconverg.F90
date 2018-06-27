!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/GenAlpha/sm_GenAlpha_checkconverg
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
!! Checks for convergence on Predictor-Corrector schemes with Generalized Alpha method
!!
!!***
#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_GenAlpha_checkconverg(testflg,ibd,convflag_pc,errmax_pc)

  use Driver_interface, only : Driver_abortFlash
  use sm_GenAlpha_data, only : sm_GenAlpha_Info, sm_GenAlpha_type

  implicit none
  
  ! IO
  integer, intent(IN)  :: testflg,ibd
  integer, intent(OUT) :: convflag_pc
  real, intent(OUT)    :: errmax_pc

  ! Internal
  type(sm_GenAlpha_type), pointer :: integ

  integ => sm_GenAlpha_info(ibd)
  convflag_pc = -99999
  errmax_pc   =  1.

  select case(testflg)
  case(SM_TESTCONVERGE) 
    convflag_pc = integ%pcconvflag
    errmax_pc   = integ%pcerr

  case(SM_SETNOTCONVERGED) ! Case one body sub-iteration has not converged.
    integ%pcconvflag = SM_PCNOTCONVERGED
    integ%pcflag     = SM_PCCORRECTOR

  case default
     call Driver_abortFlash('sm_GenAlpha_checkconverg: Unknown test flag.')
  
  end select

  return

end subroutine sm_GenAlpha_checkconverg
