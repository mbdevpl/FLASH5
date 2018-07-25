!!****f* source/Grid/GridMain/AMR/Amrex/Grid_restrictAllLevels
!!
!! NAME
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!!  call Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid-managed, cell-centered data from leaf blocks down
!!  to all refinement levels. Normally FLASH only evolves on the leaf blocks,
!!  calling this routine makes all levels have valid data.  This is mostly for
!!  visualization purposes to be able to look at different levels of resolution.
!!
!! ARGUMENTS
!!  None
!!
!!***

#include "constants.h"

subroutine Grid_restrictAllLevels()
    use Grid_data,         ONLY : gr_convertToConsvdInMeshInterp
    use gr_amrexInterface, ONLY : gr_restrictAllLevels

    implicit none

    call gr_restrictAllLevels(CENTER, &
                              convertPtoC=gr_convertToConsvdInMeshInterp, &
                              convertCtoP=gr_convertToConsvdInMeshInterp) 
end subroutine Grid_restrictAllLevels

