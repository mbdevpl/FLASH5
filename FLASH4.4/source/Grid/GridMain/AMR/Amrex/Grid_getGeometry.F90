!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getGeometry
!!
!! NAME
!!  Grid_getGeometry
!!
!! SYNOPSIS
!!
!!  Grid_getGeometry(integer (OUT)  :: geometry)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the global grid geometry.
!!  valid values are (CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL) defined
!!  in file "constants.h"
!!
!!
!! ARGUMENTS
!!
!!  geometry - returned value
!!
!!***

subroutine Grid_getGeometry(geometry)
  implicit none

#include "constants.h"

  integer, intent(OUT) :: geometry

  ! DEVNOTE: FIXME AMReX's Geometry class inherits a method called Coord() that should
  ! give this value.  Not found in Fortran interface yet.
  geometry = CARTESIAN
end subroutine Grid_getGeometry

