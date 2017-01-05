!!****if* source/diagnostics/ProtonEmission/localAPI/pem_statisticalXYZunitSphere
!!
!! NAME
!!
!!  pem_statisticalXYZunitSphere
!!
!! SYNOPSIS
!!
!!  call pem_statisticalXYZunitSphere (integer, intent (in)  :: arraySize,
!!                                     integer, intent (out) :: nSphere,
!!                                     real,    intent (out) :: xSphere (1:arraySize),
!!                                     real,    intent (out) :: ySphere (1:arraySize),
!!                                     real,    intent (out) :: zSphere (1:arraySize))
!!
!! DESCRIPTION
!!
!!  Returns a specific number of statistical (x,y,z) triples on the surface of a sphere with
!!  unit radius. The random number generator has to be inititalized (seeded) in order to be
!!  able to use this routine. The supplied (x,y,z) arrays will be used first to collect a
!!  sequence of random number triples. The (x,y,z) triples are analyzed and those inside the
!!  unit sphere are placed at the beginning of the (x,y,z) arrays. Finally, all (x,y,z) triples
!!  inside the sphere are projected (normalized) to the sphere's surface.
!!
!! ARGUMENTS
!!
!!  arraySize     : maximum number of statistical (x,y,z) triples that can be returned per call
!!  nSphere       : the actual number of statistical (x,y,z) sphere triples returned
!!  xSphere       : the x-coordinates of the statistical (x,y,z) sphere triples
!!  ySphere       : the y-coordinates of the statistical (x,y,z) sphere triples
!!  zSphere       : the z-coordinates of the statistical (x,y,z) sphere triples
!!
!! NOTES
!!
!!  The current algorithm uses the sphere (x,y,z) arrays for harvesting the random numbers needed.
!!  Since some of these numbers lay not within the sphere, the number of accepted random numbers
!!  is approximately pi/6 or 1/2 less than the maximum number of random numbers harvested.
!!  In order for this procedure to be efficient in collecting all the sphere triples, the size of
!!  the (x,y,z) arrays should not be too small.
!!
!!  Needs the (seeded) random number generator as provided in the FLASH untilities directory.
!!
!!***

subroutine pem_statisticalXYZunitSphere (arraySize,          &
                                                    nSphere, &
                                                    xSphere, &
                                                    ySphere, &
                                                    zSphere  )

  implicit none

  integer, intent (in)  :: arraySize
  integer, intent (out) :: nSphere
  real,    intent (out) :: xSphere (1:arraySize)
  real,    intent (out) :: ySphere (1:arraySize)
  real,    intent (out) :: zSphere (1:arraySize)

  nSphere = 0
  xSphere (:) = 0.0
  ySphere (:) = 0.0
  zSphere (:) = 0.0

  return
end subroutine pem_statisticalXYZunitSphere
