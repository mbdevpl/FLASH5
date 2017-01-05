!!****if* source/diagnostics/ProtonEmission/localAPI/pem_recordProtonOnScreen
!!
!! NAME
!!
!!  pem_recordProtonOnScreen
!!
!! SYNOPSIS
!!
!!  call pem_recordProtonOnScreen (real (in) :: px,
!!                                 real (in) :: py,
!!                                 real (in) :: pz,
!!                                 real (in) :: vx,
!!                                 real (in) :: vy,
!!                                 real (in) :: vz)
!!
!! DESCRIPTION
!!
!!  Records the proton on the nearest emission detector screen. The proton has an initial
!!  position and direction (velocity). The routine records the (x,y) location on
!!  the screen in terms of local screen coordinates. The (x,y) screen coordinates are
!!  rescaled to the square screen side length, such that all protons within the screen
!!  have screen coordinate paris (x,y) within the range [0,1]. Protons not hitting
!!  any screen are discarded.
!!
!!  Procedure:
!!
!!                                    |\
!!                                    | \
!!                                    |  \                   C = detector center
!!                                    |   \
!!                                    | uY|\                 n = normal unit vector
!!                                    |   | \
!!                                    |   |  \         uX,uY,n = form a right hand rule
!!                      detector      |   |   |                  orthogonal vector system:
!!                                    |   C---|---> n                  n = uY x uX
!!                       screen       \  / \  |
!!                                     \/   \ |              P = position of proton
!!                                     /\  uX\|
!!                                    /  \    |              v = velocity vector of proton
!!                                   /    \   |
!!                                  /      \  |          (x,y) = local screen coordinates
!!                                 /        \ |                  (basis = uX,uY) where proton
!!                                /          \|                   will hit the screen plane
!!                               /
!!                              /
!!                             /
!!                            P---> v -----------(x,y)
!!
!!
!!  From geometry:
!!
!!                    i) h = [(C - P) dot n] / [v dot n] > 0  ->  proton hits screen plane
!!
!!                   ii) x = uY dot [v cross (C - P)] / [v dot n]
!!
!!                  iii) y = uX dot [(C - P) cross v] / [v dot n]
!!
!!                   iv) if [v dot n] = 0  ->  parallel to screen plane, no hit
!!
!!                    v) hv -> vector between P and screen hitting position
!!
!!                   vi) min (h > 0) -> indicates nearest detector screen
!!
!! ARGUMENTS
!!
!!  px       : position global x-coordinate of the proton
!!  py       : position global y-coordinate of the proton
!!  pz       : position global z-coordinate of the proton
!!  vx       : velocity x-coordinate component of the proton
!!  vy       : velocity y-coordinate component of the proton
!!  vz       : velocity z-coordinate component of the proton
!!
!! NOTES
!!
!!***

subroutine pem_recordProtonOnScreen (px, py, pz, vx, vy, vz)

  implicit none

  real, intent (in) :: px, py, pz
  real, intent (in) :: vx, vy, vz

  return
end subroutine pem_recordProtonOnScreen
