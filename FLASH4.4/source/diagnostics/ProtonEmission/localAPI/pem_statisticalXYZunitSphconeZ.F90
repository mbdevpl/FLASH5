!!****if* source/diagnostics/ProtonEmission/localAPI/pem_statisticalXYZunitSphconeZ
!!
!! NAME
!!
!!  pem_statisticalXYZunitSphconeZ
!!
!! SYNOPSIS
!!
!!  call pem_statisticalXYZunitSphconeZ (real,    intent (in)  :: cosAlpha,
!!                                       real,    intent (in)  :: sinAlpha,
!!                                       integer, intent (in)  :: arraySize,
!!                                       integer, intent (in)  :: nWanted
!!                                       integer, intent (out) :: nReturned,
!!                                       real,    intent (out) :: xSphcone (1:arraySize),
!!                                       real,    intent (out) :: ySphcone (1:arraySize),
!!                                       real,    intent (out) :: zSphcone (1:arraySize))
!!
!! DESCRIPTION
!!
!!  Returns a specific number of statistical (x,y,z) triples laying on the surface of
!!  a unit radius spherical cone centered around the z-axis. The alignment of the unit
!!  spherical cone with the z-axis makes it easy to formulate the circumscribing rectangular
!!  prism in terms of simple x,y,z limits:
!!
!!                                      |z
!!                                      |       <- zmax = 1
!!                               \      |      /
!!                                \alpha|alpha/
!!                                 \    |    /
!!                                  \   |   /
!!                                   \  |  /
!!                                    \ | /
!!                                     \|/
!!                       -------|-------O-------|--------> x
!!                           -xmax           +xmax
!!
!!  We have the following values for the rectangular prism boundaries (a = alpha):
!!
!!                     a  |  xmax |  xmin |  ymax |  ymin |  zmax |  zmin
!!                  ------------------------------------------------------
!!                  <= 90 | sin a | -xmax | sin a | -ymax |   1   |   0
!!                   > 90 |   1   |   -1  |   1   |   -1  |   1   | cos a
!!
!!  The worst volume ratio prism/sphcone we have for a -> 0, in which case that ratio is
!!  12/pi. This means that roughly 1/4 of all statistical (x,y,z) triples will be accepted.
!!  The conversion of statistical [0,+1] numbers to statistical numbers needed in the range
!!  of the r-rescaled rectangular prism is as follows:
!!
!!     For x-coordinate: [0,+1] -> expand to [-1,+1] -> xmax * [-1,+1] -> [-xmax, xmax]
!!     For y-coordinate: [0,+1] -> expand to [-1,+1] -> ymax * [-1,+1] -> [-ymax, ymax]
!!     For z-coordinate: [0,+1] -> expand to [zmin,+1]
!!
!!  The left-end expansion from [0,+1] to [N,+1], where N is any number < 0, is achieved via:
!!
!!                   [0,+1]  ->  [0,+1] + (1 - [0,+1]) * N  =  [N,+1]
!!
!!  Hence we have:
!!                   [0,+1]  ->  [0,+1] + [0,+1] - 1           =  [-1,+1]
!!                   [0,+1]  ->  [0,+1] + (1 - [0,+1]) * zmin  =  [zmin,+1]
!!
!!  A statistical (x,y,z) triple is accepted, if the following two conditions hold:
!!
!!            Condition  I:   0 < x^2 + y^2 + z^2 <= 1
!!            Condition II:   z / sqrt (x^2 + y^2 + z^2) >= cos a
!!
!!  Condition II would involve the evaluation of a square root for every (x,y,z) triple
!!  created. To avoid this and the division operator as well, we reformulate condition II
!!  as follows:
!!
!!            Condition II:   sign (z) * z^2 >= sign (cos a) * (cos a)^2 * (x^2 + y^2 + z^2)
!!
!!  If both conditions are satisfied, the (x,y,z) triple is projected to the surface of
!!  the spherical cone. This is done by normalizing the x,y,z components.
!!
!! ARGUMENTS
!!
!!  cosAlpha      : cosine of half the spherical cone apex angle
!!  sinAlpha      : sine of half the spherical cone apex angle
!!  arraySize     : maximum number of statistical (x,y,z) triples that can be returned per call
!!  nWanted       : the wanted number of statistical (x,y,z) triples
!!  nReturned     : the actual number of statistical (x,y,z) spherical cone triples returned
!!  xSphcone      : the x-coordinates of the statistical (x,y,z) spherical cone triples
!!  ySphcone      : the y-coordinates of the statistical (x,y,z) spherical cone triples
!!  zSphcone      : the z-coordinates of the statistical (x,y,z) spherical cone triples
!!
!! NOTES
!!
!!  The current algorithm uses the sphcone (x,y,z) arrays for harvesting the random numbers needed.
!!  Since some of these numbers lay not within the spherical cone, the number of accepted random
!!  numbers is at its worst approximately pi/12 or 1/4 the maximum number of random numbers harvested.
!!  In order for this procedure to be efficient in collecting all the spherical cone triples, the
!!  size of the (x,y,z) arrays should not be too small.
!!
!!  Needs the (seeded) random number generator as provided in the FLASH untilities directory.
!!  The random number generator has to be inititalized (seeded) in order to be able to use
!!  this routine.
!!
!!  The routine does not check, if cosAlpha and sinAlpha are conmensurate. Only trivial range
!!  checks are performed. The reason for passing both these trig values in the argument is that
!!  for applications usually the cone half apex angle is fixed throughout an entire simulation.
!!  These trig values can hence be determined during initialization, avoiding repeated calls
!!  to (expensive) trig functions.
!!
!!***

subroutine pem_statisticalXYZunitSphconeZ (cosAlpha, sinAlpha,       &
                                           arraySize,                &
                                           nWanted,                  &
                                                          nReturned, &
                                                          xSphcone,  &
                                                          ySphcone,  &
                                                          zSphcone   )

  implicit none

  real,    intent (in)  :: cosAlpha, sinAlpha
  integer, intent (in)  :: arraySize
  integer, intent (in)  :: nWanted
  integer, intent (out) :: nReturned
  real,    intent (out) :: xSphcone (1:arraySize)
  real,    intent (out) :: ySphcone (1:arraySize)
  real,    intent (out) :: zSphcone (1:arraySize)

  nReturned = 0
  xSphcone (:) = 0.0
  ySphcone (:) = 0.0
  zSphcone (:) = 0.0

  return
end subroutine pem_statisticalXYZunitSphconeZ
