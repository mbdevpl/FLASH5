!!****if* source/diagnostics/ProtonEmission/localAPI/pem_activeEmissionRegion
!!
!! NAME
!!
!!  pem_activeEmissionRegion
!!
!! SYNOPSIS
!!
!!  pem_activeEmissionRegion (real (in) :: x,
!!                            real (in) :: y,
!!                            real (in) :: z)
!!
!! DESCRIPTION
!!
!!  This logical function returns true, if the (x,y,z) triple lays within at least one
!!  of the emission boxes and within the calculation domain. If no emission boxes have
!!  been specified, the entire domain is assumed to emit protons.
!!
!! ARGUMENTS
!!
!!  x : emission location x-coordinate
!!  y : emission location y-coordinate
!!  z : emission location z-coordinate
!!
!! NOTES
!!        
!!  none
!!
!!***

logical function pem_activeEmissionRegion (x,y,z)

  implicit none

  real, intent (in) :: x,y,z

  pem_activeEmissionRegion = .false.

  return
end function pem_activeEmissionRegion
