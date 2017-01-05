!!****if* source/diagnostics/ProtonEmission/localAPI/pem_protonSpeed
!!
!! NAME
!!
!!  pem_protonSpeed
!!
!! SYNOPSIS
!!
!!  pem_protonSpeed (real (in) :: EMeV)
!!
!! DESCRIPTION
!!
!!  Given the proton energy E in MeV, this function calculates the resulting
!!  relativistic proton speed:
!!
!!               v  =  c * sqrt (1.0 - (mc2 / (E + mc2)) ** 2) 
!!
!!  where E is in erg, m is the proton mass in g, c is the light speed in cm/s.
!!
!! ARGUMENTS
!!
!!  EMeV : the proton energy in MeV
!!
!! NOTES
!!        
!!  none
!!
!!***

real function pem_protonSpeed (EMeV)

  implicit none

  real, intent (in) :: EMeV

  pem_protonSpeed = 0.0

  return
end function pem_protonSpeed
