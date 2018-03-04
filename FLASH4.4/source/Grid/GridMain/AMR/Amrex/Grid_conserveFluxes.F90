!!****f* source/Grid/GridMain/AMR/Amrex/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  Grid_conserveFluxes(integer(IN) :: axis,
!!                      integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!  axis - the only acceptable value for AMReX is ALLDIR.
!!  level - refinement level whose flux data should be overwritten
!!
!!***

#include "constants.h"

subroutine Grid_conserveFluxes(axis, level)
    use Driver_interface,     ONLY : Driver_abortFlash
    use gr_physicalMultifabs, ONLY : fluxes, &
                                     flux_registers

    implicit none

    integer, intent(IN) :: axis
    integer, intent(IN) :: level

    if (axis /= ALLDIR) then
        call Driver_abortFlash("[Grid_conserveFluxes] AMReX requires axis==ALLDIR")
    end if

    ! DEV: TODO Write this once the necessary routine is available
    ! in the AMReX Fortran interface
    ! call flux_registers(level)%overwriteFluxes(fluxes(level, :))

    call Driver_abortFlash("[Grid_conserveFluxes] Not implemented yet for AMReX")
end subroutine Grid_conserveFluxes

