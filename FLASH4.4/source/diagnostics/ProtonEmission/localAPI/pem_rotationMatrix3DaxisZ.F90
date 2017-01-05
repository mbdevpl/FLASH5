!!****if* source/diagnostics/ProtonEmission/localAPI/pem_rotationMatrix3DaxisZ
!!
!! NAME
!!
!!  pem_rotationMatrix3DaxisZ
!!
!! SYNOPSIS
!!
!!  call pem_rotationMatrix3DaxisZ (real, intent (in)  :: Vx,
!!                                  real, intent (in)  :: Vy,
!!                                  real, intent (in)  :: Vz,
!!                                  real, intent (out) :: R11,
!!                                  real, intent (out) :: R12,
!!                                  real, intent (out) :: R13,
!!                                  real, intent (out) :: R21,
!!                                  real, intent (out) :: R22,
!!                                  real, intent (out) :: R23,
!!                                  real, intent (out) :: R31,
!!                                  real, intent (out) :: R32,
!!                                  real, intent (out) :: R33)
!!
!! DESCRIPTION
!!
!!  Given a vector V with 3D cartesian components (Vx,Vy,Vz), this routine establishes
!!  the 3x3 rotation matrix R, such that W = R*V results in a vector W colinear with the
!!  z axis and with equal magnitude as V, i.e |W| = |V|. The rotation matrix is calculated
!!  from 2 rotation angles theta and phi, representing respectively first a rotation around
!!  the y axis followed by a rotation around the z axes:
!!
!!                   cos(phi)cos(theta)    sin(phi)cos(theta)    - sin(theta)
!!
!!          R   =       - sin(phi)             cos(theta)              0
!!
!!                   cos(phi)sin(theta)    sin(phi)sin(theta)      cos(theta)
!!
!!  The matrix R is orthonormal and therefore R * R(transpose) = I.
!!
!! ARGUMENTS
!!
!!  Vx  : x-coordinate of vector V
!!  Vy  : y-coordinate of vector V
!!  Vz  : z-coordinate of vector V
!!  R11 : 1,1 element of rotation matrix
!!  R12 : 1,2 element of rotation matrix
!!  R13 : 1,3 element of rotation matrix
!!  R21 : 2,1 element of rotation matrix
!!  R22 : 2,2 element of rotation matrix
!!  R23 : 2,3 element of rotation matrix
!!  R31 : 3,1 element of rotation matrix
!!  R32 : 3,2 element of rotation matrix
!!  R33 : 3,3 element of rotation matrix
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pem_rotationMatrix3DaxisZ (Vx,Vy,Vz,  R11,R12,R13,R21,R22,R23,R31,R32,R33)

  implicit none

  real, intent (in)  :: Vx, Vy, Vz
  real, intent (out) :: R11,R12,R13
  real, intent (out) :: R21,R22,R23
  real, intent (out) :: R31,R32,R33

  R11 = 0.0
  R12 = 0.0
  R13 = 0.0
  R21 = 0.0
  R22 = 0.0
  R23 = 0.0
  R31 = 0.0
  R32 = 0.0
  R33 = 0.0

  return
end subroutine pem_rotationMatrix3DaxisZ
