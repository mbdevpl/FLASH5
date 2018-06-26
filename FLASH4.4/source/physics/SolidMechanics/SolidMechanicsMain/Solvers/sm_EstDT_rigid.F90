!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_EstDT_rigid
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!! 
!! Compute DT for body
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

subroutine sm_EstDT_rigid(ibd, dt)

  implicit none

  ! IO Variables
  integer, intent(in)  :: ibd
  real,    intent(out) :: dt
  
  dt = 1.e12

  return

end subroutine sm_EstDT_rigid

