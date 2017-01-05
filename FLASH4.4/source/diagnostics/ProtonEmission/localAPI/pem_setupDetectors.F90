!!****if* source/diagnostics/ProtonEmission/localAPI/pem_setupDetectors
!!
!! NAME
!!
!!  pem_setupDetectors
!!
!! SYNOPSIS
!!
!!  call pem_setupDetectors ()
!!
!! DESCRIPTION
!!
!!  Sets up the detectors (location, properties) to be used for proton emission.
!!  Since the emission screen detector files need to be open during the entire
!!  simulation, this routine opens all detector files.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All needed info is read in as runtime parameters.
!!
!!***

subroutine pem_setupDetectors ()

  implicit none

  return
end subroutine pem_setupDetectors
