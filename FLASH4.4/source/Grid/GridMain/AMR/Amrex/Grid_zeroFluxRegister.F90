!!****if* source/Grid/AMR/Amrex/Grid_zeroFluxRegister
!!
!! NAME
!!  Grid_zeroFluxRegister
!!
!! SYNOPSIS
!!  call Grid_zeroFluxRegister(integer(IN) :: fine_level)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  Given an
!!  index for the fine level, set the fine and coarse data to zero in the
!!  associated flux register.
!!
!! ARGUMENTS
!!  fine_level - the 1-based level index (1 is the coarsest level) that 
!!               indicates the fine level of the flux register on which to 
!!               operate.
!!
!! SEE ALSO
!!   Grid_addFineToFluxRegister
!!   Grid_addCoarseToFluxRegister
!!
!!***

#include "Flash.h"

subroutine Grid_zeroFluxRegister(fine_level)
    use amrex_fort_module,    ONLY : wp => amrex_real
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    use Driver_interface,     ONLY : Driver_abortFlash
    use gr_physicalMultifabs, ONLY : flux_registers

    implicit none

    integer, intent(IN) :: fine_level

    integer :: fine

    ! Skip if we aren't using fluxes
    if (NFLUXES < 1) then
        RETURN
    end if

    ! FLASH uses 1-based level index / AMReX uses 0-based index
    fine = fine_level - 1

    if ((fine <= 0) .OR. (fine > amrex_get_finest_level())) then
        call Driver_abortFlash("[Grid_zeroFluxRegister] Invalid level")
    end if

    call flux_registers(fine)%setval(0.0_wp)
end subroutine Grid_zeroFluxRegister

