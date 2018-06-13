!!****if* source/Grid/GridMain/AMR/Amrex/Grid_putFluxData
!!
!! NAME
!!  Grid_putFluxData
!!
!! SYNOPSIS
!!
!!  call Grid_putFluxData(integer(IN)         :: level,
!!                        integer(IN)         :: axis,
!!              optional, integer(IN), target :: pressureSlots(:),
!!              optional, real(IN),           :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!  Request that the Grid unit load into its flux registers the flux data
!!  corresponding to the given axial directions and defined on the given
!!  refinement level.  This routine needs to be used with adaptive mesh
!!  since fluxes calculated by two blocks at a fine/coarse boundary have 
!!  different accuracy.
!!
!!  It is assumed that before calling this routine, the client code has already
!!  written flux data to Grid's data structures using the Grid_getFluxPtr
!!  interface.  Once data is loaded in the registers, a call to
!!  Grid_conserveFluxes applies the flux conservation algorithm to make it
!!  consistent across the fine/coarse boundaries.
!!
!! ARGUMENTS
!!  level - the 1-based level index (1 is the coarsest level) indicating which
!!          refinement level's fluxes should be loaded in the flux registers.
!!  axis - integer value specifying on which cell faces to put fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!  pressureSlots - If present and greater than zero, this indicates one or more flux
!!                  variables in the fluxes array that may need special handling because
!!                  they really scale like flux densities when other flux variables scale
!!                  like fluxes.  For the split PPM Hydro implementation, for example,
!!                  this is normally used for the pressure "flux" variable, but
!!                  it can be applied to other flux variables that the caller keeps in
!!                  flux density form.
!!                  The special handling consists in multiplying the flux variables with
!!                  the appropriate face areas, which are taken from the areaLeft array
!!                  argument, before storing them in the arrays on which the Grid unit
!!                  acts.
!!                  Special handling should only be requested if the Grid units handles
!!                  flux variable "as fluxes" (or else if it does nto matter anyway,
!!                  as is the case in 1D).
!!                  The pressureSlots argument in the corresponding Grid_getFluxData
!!                  call should generally match the one in the Grid_putFluxData call.
!!  areaLeft - areas of left and right faces, only used if special scaling is
!!             requested with the pressureSlot argument.
!!
!! SEE ALSO
!!   Grid_getFluxPtr/Grid_releaseFluxPtr
!!   Grid_conserveFluxes
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_putFluxData(level, axis, pressureSlots, areaLeft)
  use amrex_fort_module,    ONLY : wp => amrex_real

  use Driver_interface,     ONLY : Driver_abortFlash
  use Grid_interface,       ONLY : Grid_getGeometry
  use gr_physicalMultifabs, ONLY : flux_registers, &
                                   fluxes

  implicit none

  integer,                intent(IN)                   :: level
  integer,                intent(IN), optional         :: axis
  integer,                intent(IN), optional, target :: pressureSlots(:)
  real,                   intent(IN), optional         :: areaLeft(:,:,:)

  integer :: geometry

  if (present(axis)) then
    call Driver_abortFlash("[Grid_putFluxData] axis not accepted with AMReX")
  end if
  if (present(pressureSlots)) then
    call Driver_abortFlash("[Grid_putFluxData] pressureSlots not accepted with AMReX")
  end if
  if (present(areaLeft)) then
    call Driver_abortFlash("[Grid_putFluxData] areaLeft not accepted with AMReX")
  end if

  if(NFLUXES < 1)   RETURN

  ! FLASH uses 1-based level index / AMReX uses 0-based index
  call flux_registers(level-1)%setval(0.0_wp)
  
  call Grid_getGeometry(geometry)
  
  select case (geometry)
  case (CARTESIAN)
    ! DEV: TODO This routine should take a densityMask array as an argument
    ! so that the routine can determine which need to be scaled and which do not
#if   NDIM == 2
    call flux_registers(level-1)%fineadd(fluxes(level-1, 1:NDIM), 0.5_wp)
#elif NDIM == 3
    call flux_registers(level-1)%fineadd(fluxes(level-1, 1:NDIM), 0.25_wp)
#endif
  case default
    call Driver_abortFlash("[Grid_putFluxData] Only works with Cartesian")
  end select
end subroutine Grid_putFluxData

