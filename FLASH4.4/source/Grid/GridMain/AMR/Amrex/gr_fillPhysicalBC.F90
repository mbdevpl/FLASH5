!!****if* source/Grid/GridMain/AMR/Amrex/gr_fillPhysicalBC
!!
!! NAME
!!
!!  gr_fillPhysicalBC
!!
!! SYNOPSIS
!!
!!  call gr_fillPhysicalBC(amrex_multifab(IN) :: pmf,
!!                         integer(IN)        :: scomp,
!!                         integer(IN)        :: ncomp,
!!                         amrex_real(IN)     :: time,
!!                         amrex_geometry(IN) :: pgeom)
!!
!! DESCRIPTION 
!!  
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
!!
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

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine gr_fillPhysicalBC(pmf, scomp, ncomp, time, pgeom) bind(c)
    use iso_c_binding
    
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom
    use amrex_box_module,       ONLY : amrex_box
    use amrex_geometry_module,  ONLY : amrex_geometry, &
                                       amrex_is_all_periodic
    use amrex_multifab_module,  ONLY : amrex_multifab, &
                                       amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy

    use Driver_interface,       ONLY : Driver_abortFlash
    use Grid_data,              ONLY : gr_maxRefine, &
                                       gr_domainBC
    use Grid_interface,         ONLY : Grid_bcApplyToRegion, &
                                       Grid_bcApplyToRegionSpecialized
    use gr_amrexInterface,      ONLY : gr_splitFabAtBoundary, &
                                       gr_copyFabInteriorToRegion, &
                                       gr_copyGuardcellRegionToFab
    use block_metadata,         ONLY : block_metadata_t 

    implicit none

    type(c_ptr),    value :: pmf
    type(c_ptr),    value :: pgeom
    integer(c_int), value :: scomp
    integer(c_int), value :: ncomp
    real(wp),       value :: time

    type(amrex_geometry)   :: geom
    type(amrex_multifab)   :: mfab
    type(amrex_mfiter)     :: mfi
    type(amrex_box)        :: box
    type(block_metadata_t) :: blockDesc

    integer :: j
    real    :: delta, delta_j
    integer :: level, face, dir
    integer :: finest_level
    integer :: axis, axis2, axis3

    logical :: mask(ncomp)
    logical :: applied

    integer                       :: endPts(LOW:HIGH, 1:MDIM)
    integer                       :: regionType(1:MDIM)
    integer                       :: regionSize(4)
    real(wp), pointer, contiguous :: regionData(:, :, :, :)
    
    real(wp), pointer, contiguous :: solnData(:, :, :, :)
    integer                       :: limitsGC(LOW:HIGH, 1:MDIM)
    integer                       :: interior(LOW:HIGH, 1:MDIM)
    integer                       :: guardcells(LOW:HIGH, 1:MDIM)

    logical :: found

    geom = pgeom
    mfab = pmf

    ! When fillpatch is called, it is passed data for a known level.
    ! However, this means that fillpatch does not have knowledge of the
    ! level number.  Therefore, when fillpatch calls this routine, we as well
    ! do not know the level number.
    !
    ! DEV: FIXME User BC routines could need the level information to obtain the
    ! deltas and to determine which multifab to get data from.  We have at least
    ! the following options:
    !   1) Keep this hack
    !   2) Get AMReX to give level value instead of pgeom
    !   3) Get AMReX to include level in amrex_geometry
    delta = geom%dx(IAXIS)
    level = INVALID_LEVEL

    ! AMReX uses 0-based level index set / FLASH uses 1-based
    do j = 0, gr_maxRefine
        delta_j = amrex_geom(j)%dx(IAXIS)
        if (delta == delta_j) then
            level = j + 1
            EXIT
        end if
    end do

    if (level == INVALID_LEVEL) then
        call Driver_abortFlash("[gr_fillPhysicalBC] Could not reverse engineer level")
    end if

#ifdef DEBUG_GRID
    write(*,'(A,A,I3)') "[gr_fillPhysicalBC]", &
                        "                  Level ", level
#endif

    call amrex_mfiter_build(mfi, mfab, tiling=.false.)
    do while(mfi%next())
       ! 0-based, cell-centered, and global indices
       solnData => mfab%dataPtr(mfi)

       ! 0-based, cell-centered, and global indices
       box = mfi%fabbox()
       limitsGC(:, :) = 0
       limitsGC(LOW,  1:NDIM) = box%lo(1:NDIM)
       limitsGC(HIGH, 1:NDIM) = box%hi(1:NDIM)

       ! DEV: The given box is not necessarily a FLASH block, but rather
       ! could be an arbitrary rectangular region in the domain and its GC
       !
       ! One result of this is that grid_index is non-sensical here
       blockDesc%level = level
       blockDesc%grid_index = -1
       blockDesc%limits(LOW,  :)   = limitsGC(LOW,  :) + 1 + NGUARD
       blockDesc%limits(HIGH, :)   = limitsGC(HIGH, :) + 1 - NGUARD
       blockDesc%limitsGC(LOW,  :) = limitsGC(LOW,  :) + 1
       blockDesc%limitsGC(HIGH, :) = limitsGC(HIGH, :) + 1

       ! Check for boundaries on both faces along all directions
       do axis = 1, NDIM
          ! Enforce index ordering needed for region/Grid_bcApplyToRegion
          if      (axis == IAXIS) then
             axis2 = JAXIS
             axis3 = KAXIS
          else if (axis == JAXIS) then
             axis2 = IAXIS
             axis3 = KAXIS
          else
             axis2 = IAXIS
             axis3 = JAXIS
          end if

          do face = LOW, HIGH
             ! We have configured AMReX to handle periodic BC automatically
             ! DEV FIXME: This does not allow users to provide custom periodic
             ! BC code.
             if (gr_domainBC(face, axis) == PERIODIC)    CYCLE

             ! interior/guardcells are 0-based, cell-centered, and global indices
             call gr_splitFabAtBoundary(face, axis, limitsGC, &
                                        amrex_geom(level-1)%dx, &
                                        interior, guardcells, found)
             if (.NOT. found)   CYCLE

             ! Grow GC region to meet needs of Grid_bcApplyToRegion
             ! i.e. include NGUARD cells on either side of boundary and 1-based
             endPts(:, :) = 1
             endPts(:, 1:NDIM) = guardcells(:, 1:NDIM) + 1
             if (face == LOW) then
                endPts(LOW,  axis) = (guardcells(HIGH, axis) + 1) - (NGUARD - 1)
                endPts(HIGH, axis) = (guardcells(HIGH, axis) + 1) +  NGUARD
             else
                endPts(LOW,  axis) = (interior(HIGH, axis)   + 1) - (NGUARD - 1)
                endPts(HIGH, axis) = (interior(HIGH, axis)   + 1) +  NGUARD
             end if

             ! Create buffer to hold data with indices permuted for use with
             ! Grid_bcApplyToRegion
             regionSize(BC_DIR)     = endPts(HIGH, axis)  - endPts(LOW, axis)  + 1
             regionSize(SECOND_DIR) = endPts(HIGH, axis2) - endPts(LOW, axis2) + 1
             regionSize(THIRD_DIR)  = endPts(HIGH, axis3) - endPts(LOW, axis3) + 1
             regionSize(STRUCTSIZE) = ncomp

             ! 1-based, cell-centered, and local indices 
             allocate(regionData(regionSize(BC_DIR), &
                                 regionSize(SECOND_DIR), &
                                 regionSize(THIRD_DIR), &
                                 regionSize(STRUCTSIZE)) )

             regionData(:, :, :, :) = 0.0d0
             call gr_copyFabInteriorToRegion(solnData, face, axis, &
                                             interior, scomp, ncomp, regionData)

             ! As regionData only contains those physical quantities that AMReX
             ! asks for, no need for masking
             mask = .TRUE.

             ! Let simulation do BC fill if so desired
             applied = .FALSE.
             call Grid_bcApplyToRegionSpecialized(gr_domainBC(face, axis), &
                                                  CENTER, NGUARD, &
                                                  axis, face, &
                                                  regionData, regionSize, &
                                                  mask, applied, blockDesc, &
                                                  axis2, axis3, endPts, 0)

             if (.NOT. applied) then
                ! Have FLASH fill GC in special data buffer
                call Grid_bcApplyToRegion(gr_domainBC(face, axis), &
                                          CENTER, NGUARD, &
                                          axis, face, &
                                          regionData, regionSize, &
                                          mask, applied, blockDesc, &
                                          axis2, axis3, endPts, 0)
             end if

             call gr_copyGuardcellRegionToFab(regionData, face, axis, &
                                              guardcells, scomp, ncomp, solnData)

             deallocate(regionData)

             if (.NOT. applied) then
                call Driver_abortFlash("[gr_fillPhysicalBC] BC not applied")
             end if
          end do
       end do

       nullify(solnData)
    end do
    
    call amrex_mfiter_destroy(mfi)
end subroutine gr_fillPhysicalBC

