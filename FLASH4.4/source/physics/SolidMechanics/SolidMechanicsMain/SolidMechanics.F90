!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SolidMechanics
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
!!
!!
!!
!!***
#include "SolidMechanics.h"

subroutine SolidMechanics(selector_flag, restart, convflag_all)
  use Driver_interface,  only: Driver_abortFlash
  use sm_integinterface, only: sm_integ_advance, sm_integ_adv1dt, &
                               sm_integ_checkconverg
  implicit none

  ! IO Vars
  integer, intent(in)              :: selector_flag
  logical, optional, intent(in)    :: restart
  integer, optional, intent(inout) :: convflag_all 

  ! Internal Variables
  logical :: restart_local = .false.
 
  if( present(restart) ) restart_local = restart

  select case ( selector_flag )

  case( SM_ADVANCE )

     ! Call wrapper routine for advance
     call sm_integ_advance(restart_local)

  case( SM_CHECKCONVERG )

     ! Call wrapper to check all bocies are converged 
     call sm_integ_checkconverg(convflag_all) ! This datum is also stored in sm_convflag_all
                                              ! in sm_integdata.F90

  case( SM_ADVANCE1DT )

     ! Call wrapper to iterate to convergence the Bodies
     ! 1 timestep
     call sm_integ_adv1dt(restart_local)

  case( SM_WRITECHECKPT )

     ! wrapper for writing out the checkpoint files
     call sm_integ_writeCheckpoint()

  case default
     call Driver_abortFlash("Called SolidMechanics with an unknown input.")

  end select

  return

end subroutine SolidMechanics
