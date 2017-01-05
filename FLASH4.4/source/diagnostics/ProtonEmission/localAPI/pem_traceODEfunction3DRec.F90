!!****if* source/diagnostics/ProtonEmission/localAPI/pem_traceODEfunction3DRec
!!
!! NAME
!!
!!  pem_traceODEfunction3DRec
!!
!! SYNOPSIS
!!
!!  pem_traceODEfunction3DRec (real, intent (in) :: t,
!!                             real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the proton trace ODE function in 3D rectangular coordinates to be used for Runge Kutta
!!  proton tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (px,py,pz,vx,vy,vz   position and velocity)
!!
!!***

function pem_traceODEfunction3DRec (t,y)

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: pem_traceODEfunction3DRec (1:size (y))     ! (compatible to Runge Kutta interfaces)

  pem_traceODEfunction3DRec (1:size (y)) = 0.0

  return
end function pem_traceODEfunction3DRec
