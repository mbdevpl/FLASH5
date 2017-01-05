!!****if* source/diagnostics/ProtonEmission/localAPI/pem_detectorsCheck3DRec
!!
!! NAME
!!
!!  pem_detectorsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call pem_detectorsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected detectors data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!         1) Check, if all detector screen areas are completely outside of the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pem_detectorsCheck3DRec ()

  implicit none

  return
end subroutine pem_detectorsCheck3DRec
