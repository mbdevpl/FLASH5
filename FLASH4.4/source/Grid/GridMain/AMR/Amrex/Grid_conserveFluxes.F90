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
    use Grid_interface,       ONLY : Grid_getGeometry
    use gr_physicalMultifabs, ONLY : fluxes, &
                                     flux_registers

    implicit none

    integer, intent(IN) :: axis
    integer, intent(IN) :: level

    integer :: geometry

    if (axis /= ALLDIR) then
        call Driver_abortFlash("[Grid_conserveFluxes] AMReX requires axis==ALLDIR")
    end if

    ! The AMReX flux registers are dealing with fluxes and *not* flux densities
    call Grid_getGeometry(geometry)
    
    select case (geometry)
    case (CARTESIAN)
        ! DEV: TODO This routine should take a densityMask array as an argument
        ! so that the routine can determine which need to be scaled and which do not
        call flux_registers(level)%overwrite(fluxes(level-1, :), 1.0_wp)
    case default
        call Driver_abortFlash("[Grid_conserveFluxes] Only works with Cartesian")
    end select
end subroutine Grid_conserveFluxes

