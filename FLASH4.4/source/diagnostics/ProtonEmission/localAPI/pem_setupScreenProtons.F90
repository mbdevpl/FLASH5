!!****if* source/diagnostics/ProtonEmission/localAPI/pem_setupScreenProtons
!!
!! NAME
!!
!!  pem_setupScreenProtons
!!
!! SYNOPSIS
!!
!!  call pem_setupScreenProtons ()
!!
!! DESCRIPTION
!!
!!  Sets up the emission screen protons, which means to allocate the needed screen proton
!!  array. The screen protons are used to accumulate emission protons on the detector
!!  screen(s). The routine also sets up and commits a mpi type structure corresponding to
!!  the screen protons, which will be used when sending screen protons between processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pem_setupScreenProtons ()

  implicit none

  return
end subroutine pem_setupScreenProtons
