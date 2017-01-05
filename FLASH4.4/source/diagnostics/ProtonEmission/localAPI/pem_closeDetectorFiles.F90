!!****if* source/diagnostics/ProtonEmission/localAPI/pem_closeDetectorFiles
!!
!! NAME
!!
!!  pem_closeDetectorFiles
!!
!! SYNOPSIS
!!
!!  call pem_closeDetectorFiles ()
!!
!! DESCRIPTION
!!
!!  Closes all detector files that are currently open. Only the master processor can close
!!  detector files.
!!
!! ARGUMENTS
!!
!!***

subroutine pem_closeDetectorFiles ()

  implicit none

  return
end subroutine pem_closeDetectorFiles
