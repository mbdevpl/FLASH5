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
    logical                :: ntype(MDIM)
    integer                :: gds
    type(amrex_mfiter)     :: mfi
    type(amrex_box)        :: box
    type(block_metadata_t) :: blockDesc
    
    integer :: n_cells_domain
    integer :: n_cells_level

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
    ! Therefore, we reverse engineer the level by looking at the number of cells
    ! in the physical domain at the refinement level associated with the given
    ! multifab
    level = INVALID_LEVEL
    n_cells_domain = geom%domain%hi(IAXIS) - geom%domain%lo(IAXIS)
    do j = 0, gr_maxRefine
        n_cells_level= amrex_geom(j)%domain%hi(IAXIS) - &
                       amrex_geom(j)%domain%lo(IAXIS)
        if (n_cells_domain == n_cells_level) then
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

    ntype = mfab%nodal_type()
    if      (      (.NOT. ntype(IAXIS))  &
             .AND. (.NOT. ntype(JAXIS)) &
             .AND. (.NOT. ntype(KAXIS))) then
        gds = CENTER
    else if (             ntype(IAXIS)  &
             .AND. (.NOT. ntype(JAXIS)) &
             .AND. (.NOT. ntype(KAXIS))) then
        gds = FACEX
        write(*,*) "I am a facevarx multifab"
    else if (      (.NOT. ntype(IAXIS)) &
             .AND.        ntype(JAXIS) &
             .AND. (.NOT. ntype(KAXIS))) then
        gds = FACEY
        write(*,*) "I am a facevary multifab"
    else if (      (.NOT. ntype(IAXIS))  &
             .AND. (.NOT. ntype(JAXIS)) &
             .AND.        ntype(KAXIS)) then
        gds = FACEZ
        write(*,*) "I am a facevarz multifab"
    else
        call Driver_abortFlash("[gr_fillPhysicalFaceBC] I am some exotic " // &
                               "beast that exists outside the realm of FLASH")
    end if

    call amrex_mfiter_build(mfi, mfab, tiling=.false.)
    do while(mfi%next())
       ! 0-based and global indices
       solnData => mfab%dataPtr(mfi)

       ! 0-based and global indices
       box = mfi%fabbox()
       limitsGC(:, :) = 0
       limitsGC(LOW,  1:NDIM) = box%lo(1:NDIM)
       limitsGC(HIGH, 1:NDIM) = box%hi(1:NDIM)

       ! In FLASH, limits/limitsGC in a block descriptor must be 
       ! cell-centered
       if      (gds == FACEX) then
          limitsGC(HIGH, IAXIS) = limitsGC(HIGH, IAXIS) - 1
       else if (gds == FACEY) then
          limitsGC(HIGH, JAXIS) = limitsGC(HIGH, JAXIS) - 1
       else if (gds == FACEZ) then
          limitsGC(HIGH, KAXIS) = limitsGC(HIGH, KAXIS) - 1
       end if

       ! The given box is not necessarily a FLASH block, but rather
       ! could be an arbitrary rectangular region in the domain and its GC
       !
       ! One result of this is that grid_index is non-sensical here
       !
       ! 1-based, cell-centered and global indices
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

             ! If necessary, convert interior, guardcells, and endPts
             ! to match index space of given multifab
             if      (gds == FACEX) then
                endPts(    HIGH, IAXIS) = endPts(    HIGH, IAXIS) + 1
                interior(  HIGH, IAXIS) = interior(  HIGH, IAXIS) + 1
                guardcells(HIGH, IAXIS) = guardcells(HIGH, IAXIS) + 1
             else if (gds == FACEY) then
                endPts(    HIGH, JAXIS) = endPts(    HIGH, JAXIS) + 1
                interior(  HIGH, JAXIS) = interior(  HIGH, JAXIS) + 1
                guardcells(HIGH, JAXIS) = guardcells(HIGH, JAXIS) + 1
             else if (gds == FACEZ) then
                endPts(    HIGH, KAXIS) = endPts(    HIGH, KAXIS) + 1
                interior(  HIGH, KAXIS) = interior(  HIGH, KAXIS) + 1
                guardcells(HIGH, KAXIS) = guardcells(HIGH, KAXIS) + 1
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
             call gr_copyFabInteriorToRegion(solnData, gds, face, axis, &
                                             interior, scomp, ncomp, regionData)

             ! As regionData only contains those physical quantities that AMReX
             ! asks for, no need for masking
             mask = .TRUE.

             ! Let simulation do BC fill if so desired
             applied = .FALSE.
             call Grid_bcApplyToRegionSpecialized(gr_domainBC(face, axis), &
                                                  gds, NGUARD, &
                                                  axis, face, &
                                                  regionData, regionSize, &
                                                  mask, applied, blockDesc, &
                                                  axis2, axis3, endPts, 0)

             if (.NOT. applied) then
                ! Have FLASH fill GC in special data buffer
                call Grid_bcApplyToRegion(gr_domainBC(face, axis), &
                                          gds, NGUARD, &
                                          axis, face, &
                                          regionData, regionSize, &
                                          mask, applied, blockDesc, &
                                          axis2, axis3, endPts, 0)
             end if

             call gr_copyGuardcellRegionToFab(regionData, gds, face, axis, &
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

