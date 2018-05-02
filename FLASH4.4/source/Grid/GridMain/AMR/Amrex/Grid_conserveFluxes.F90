!!****f* source/Grid/GridMain/AMR/Amrex/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  call Grid_conserveFluxes(integer(IN) :: axis,
!!                           integer(IN) :: level)
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
    use amrex_fort_module,    ONLY : wp => amrex_real

    use Driver_interface,     ONLY : Driver_abortFlash
    use gr_physicalMultifabs, ONLY : fluxes, &
                                     flux_registers

    implicit none

    integer, intent(IN) :: axis
    integer, intent(IN) :: level

    if (axis /= ALLDIR) then
        call Driver_abortFlash("[Grid_conserveFluxes] AMReX requires axis==ALLDIR")
    end if

    ! The fluxes data structure should contain fluxes and *not* flux density
    ! Therefore, we do not need AMReX to scale for us.
    call flux_registers(level+1)%overwrite(fluxes(level, :), 1.0_wp)
end subroutine Grid_conserveFluxes

