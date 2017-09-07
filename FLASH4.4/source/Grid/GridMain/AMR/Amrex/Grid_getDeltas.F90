!!****if* source/Grid/GridMain/paramesh/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  Grid_getDeltas(integer(IN) :: blockId,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the grid spacing dx/dy/dz for a given blockId on the Grid.
!!  dx is the size of one cell in the x direction of a block.
!!  
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!***

subroutine Grid_getDeltas(lev, del)
  use amrex_amrcore_module, ONLY : amrex_geom
  
  implicit none

#include "constants.h"
  
  integer, intent(IN)  :: lev
  real,    intent(OUT) :: del(MDIM)

  ! AMReX uses zero-based level indexing, but FLASH assumes one-based
  del = 0.0d0
  del(1:MDIM) = amrex_geom(lev-1)%dx(1:MDIM)
end subroutine Grid_getDeltas

