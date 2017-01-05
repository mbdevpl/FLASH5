!!****if* source/diagnostics/ProtonEmission/localAPI/pem_statisticalSetSeed
!!
!! NAME
!!
!!  pem_statisticalSetSeed
!!
!! SYNOPSIS
!!
!!  call pem_statisticalSetSeed (integer, intent (in) :: seed)
!!
!! DESCRIPTION
!!
!!  Initializes the random number generator with the supplied seed. After calling this
!!  routine, the random number generator is ready for use.
!!
!! ARGUMENTS
!!
!!  seed : the seed value for the random number generator
!!
!!***

subroutine pem_statisticalSetSeed (seed)

  implicit none

  integer, intent (in) :: seed

  return
end subroutine pem_statisticalSetSeed
