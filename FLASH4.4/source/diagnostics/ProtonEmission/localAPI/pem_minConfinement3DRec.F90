!!****if* source/diagnostics/ProtonEmission/localAPI/pem_minConfinement3DRec
!!
!! NAME
!!
!!  pem_minConfinement3DRec
!!
!! SYNOPSIS
!!
!!  pem_minConfinement3DRec (integer, intent (in) :: nc,
!!                           real,    intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Places the minimum confinement rules for 3D rectangular ray tracing Runge Kutta runs
!!  as functions of the dependent variables array.
!!
!! ARGUMENTS
!!
!!  nc : number of confined dependent variables (must be 3 and <= size (y), not checked here!)
!!  y  : dependent variable array (not used here)
!!
!!***

function pem_minConfinement3DRec (nc,y)

  implicit none

  integer, intent (in) :: nc               ! it is absolutely mandatory to
  real,    intent (in) :: y (:)            ! declare the variables and
                                           ! array function in the way shown.
  real :: pem_minConfinement3DRec (1:nc)   ! (compatible to Runge Kutta interfaces)

  pem_minConfinement3DRec (1:nc) = 0.0

  return
end function pem_minConfinement3DRec
