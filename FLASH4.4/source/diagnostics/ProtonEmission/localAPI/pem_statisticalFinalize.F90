!!****if* source/diagnostics/ProtonEmission/localAPI/pem_statisticalFinalize
!!
!! NAME
!!
!!  pem_statisticalFinalize
!!
!! SYNOPSIS
!!
!!  call pem_statisticalFinalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the statistical environment that was used by the proton emission code. It currently
!!  consists in deallocation of the seed array that was needed for the random number generator.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pem_statisticalFinalize ()

  implicit none

  return
end subroutine pem_statisticalFinalize
