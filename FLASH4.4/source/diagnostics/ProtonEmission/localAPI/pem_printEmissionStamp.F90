!!****if* source/diagnostics/ProtonEmission/localAPI/pem_printEmissionStamp
!!
!! NAME
!!
!!  pem_printEmissionStamp
!!
!! SYNOPSIS
!!
!!  call pem_printEmissionStamp (real, intent (in) :: timeStep,
!!                               real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Prints out the current emission statistics of the proton emission unit:
!!
!!        i) total number of protons emitted so far
!!       ii) number of protons emitted during the current timestep
!!
!! ARGUMENTS
!!
!!  timeStep       : Current time step value
!!  timeSimulation : Current simulation time
!!
!! NOTES
!!
!!  Only the master processor writes to the emission profile file.
!!
!!***

subroutine pem_printEmissionStamp (timeStep, timeSimulation)

  implicit none
   
  real, intent (in) :: timeStep
  real, intent (in) :: timeSimulation

  return
end subroutine pem_printEmissionStamp
