!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_EstDT_rbc
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

subroutine sm_EstDT_rbc(ibd, dt)

  implicit none

  ! IO Variables
  integer, intent(in)  :: ibd
  real,    intent(out) :: dt

  ! At the time being the time step is fixed to the initial dt 
  dt = 1e5
  write(*,*) '*** Hack *** : change DT for rbc body'

end subroutine sm_EstDT_rbc

