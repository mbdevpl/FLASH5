!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_PredCorr_checkconverg
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
!! Stub Function.
!!***
subroutine sm_PredCorr_checkconverg(testflg,ibd,convflag_pc,errmax_pc)
  use Driver_interface, only : Driver_abortFlash
  implicit none
  integer, intent(IN)  :: testflg,ibd
  integer, intent(OUT) :: convflag_pc
  real, intent(OUT)    :: errmax_pc


  convflag_pc = -999999
  errmax_pc   = 1.e12

  call Driver_abortFlash('sm_pc_checkconverg Error: Include PredCorr directory for PC integration.')

  return

end subroutine sm_PredCorr_checkconverg
