!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!  Grid_getDeltas(integer(IN) :: level,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the size (dX, dY, dZ) of cells in the given refinement level where dX
!!  is the width of the cell along the first dimension (IAXIS).
!! 
!!  If a dimension is Cartesian, then the returned size is the length of any
!!  cell along that direction.  If a dimension is angular, then the size is the
!!  angle expressed in radians (as opposed to the arclength).
!!
!! ARGUMENTS 
!!  level - local block number
!!  del - array of size MDIM returned holding the dX, dY, and dZ values
!!
!!***

subroutine Grid_getDeltas(lev, del)
  use amrex_amrcore_module, ONLY : amrex_geom
  
  implicit none

#include "constants.h"
  
  integer, intent(IN)  :: lev
  real,    intent(OUT) :: del(MDIM)

  ! AMReX uses zero-based level indexing, but FLASH assumes one-based
  del = 0.0
  del(1:MDIM) = amrex_geom(lev-1)%dx(1:MDIM)
end subroutine Grid_getDeltas

