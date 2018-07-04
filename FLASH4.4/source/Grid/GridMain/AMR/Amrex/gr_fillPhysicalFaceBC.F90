!!****if* source/Grid/GridMain/AMR/Amrex/gr_fillPhysicalFaceBC
!!
!! NAME
!!  gr_fillPhysicalFaceBC
!!
!! SYNOPSIS
!!  call gr_fillPhysicalFaceBC(amrex_multifab(IN) :: pmf,
!!                             integer(IN)        :: scomp,
!!                             integer(IN)        :: ncomp,
!!                             amrex_real(IN)     :: time,
!!                             amrex_geometry(IN) :: pgeom)
!!
!! DESCRIPTION 
!!  This routine is a callback function that is given to AMReX when using the
!!  fillpatch routines.  It is given a multifab where each FAB already contains
!!  valid interior data and needs to have its guardcells filled with data that
!!  satisfies the boundary conditions of the problem.
!!
!!  This routine executes the GC fill using the GridBoundaryConditions subunit.
!!  In particular, client code is first given the opportunity to execute the
!!  fill via the routine Grid_bcApplyToRegionSpecialized.  If this routine does
!!  not handle the fill, then the fill is done via Grid_bcApplyToRegion.
!!
!! ARGUMENTS 
!!  pmf - the multifab on which to operate
!!  scomp - the 1-based index of the first physical quantity on which to carry
!!          out the operation
!!  ncomp - the number of physical quantities on which to carry out the
!!          operation
!!  time - not used in the FLASH implementation
!!  pgeom - an instance of amrex_geometry that indicates the level, deltas
!!          for all FABS in pmf
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_fillPhysicalFaceBC(pmf, scomp, ncomp, time, pgeom) bind(c)
    use iso_c_binding
    
    use amrex_fort_module,      ONLY : wp => amrex_real

    use Driver_interface,       ONLY : Driver_abortFlash
    use Grid_data,              ONLY : gr_domainBC

    implicit none

    type(c_ptr),    value :: pmf
    type(c_ptr),    value :: pgeom
    integer(c_int), value :: scomp
    integer(c_int), value :: ncomp
    real(wp),       value :: time

    integer :: axis
    integer :: face

    ! DEV: Not yet clear that we need a different BC routine for the face
    ! variables.  Concerns are:
    ! 1) will routines that create special data structure for GridBC subunit
    !    work equally well with facevar[xyz]
    ! 2) Does the normal applyBc routine work with facevar[xyz]?  Do we need
    !    to call it at all for facevars?
    ! 3) Is it easy/reasonable to determine if given mfab is facevar?

    ! DEV: TODO Allow for non-periodic BC to be used in conjunction with
    ! face vars.  At the very least we need to allow simulation code to handle
    ! these with the specialized routine
    do axis = 1, NDIM
        do face = LOW, HIGH
            if (gr_domainBC(face, axis) /= PERIODIC) then
                call Driver_abortFlash('[gr_fillPhysicalFaceBC] BC must be periodic with face vars')
            end if
        end do
    end do

    ! NO-OP.  AMReX handles periodic BC
end subroutine gr_fillPhysicalFaceBC

