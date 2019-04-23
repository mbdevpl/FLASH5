!!****if* source/Grid/GridMain/AMR/Amrex/gr_fillPhysicalBC
!!
!! NAME
!!  gr_fillPhysicalBC
!!
!! SYNOPSIS
!!  call gr_fillPhysicalBC(amrex_multifab(IN) :: pmf,
!!                         integer(IN)        :: scomp,
!!                         integer(IN)        :: ncomp,
!!                         amrex_real(IN)     :: time,
!!                         amrex_geometry(IN) :: pgeom)
!!
!! DESCRIPTION 
!!  This routine is a callback function that is given to AMReX when using the
!!  fillpatch routines.  It is given a multifab where each FAB already contains
!!  data that is valid in the domain and on the domain boundary for the case of
!!  face-centered variables.  If a LOW/HIGH pair of faces have periodic BC, then
!!  the layers of data outside these faces have also been filled with correct
!!  data.  It is expected that this routine will fill those guardcells outside
!!  the domain with data that satisfies the boundary conditions of the problem.
!!
!!  To be meaningful, it is assumed that all fabs that span a domain boundary
!!  have at least NGUARD layers of interior data along the axis of the
!!  boundary.  It is also assumed that these fabs do not require filling more
!!  than NGUARD layers of guardcells outside the domain boundary.
!!
!!  This routine executes the GC fill using the GridBoundaryConditions subunit.
!!  In particular, client code is first given the opportunity to execute the
!!  fill via the routine Grid_bcApplyToRegionSpecialized.  If this routine does
!!  not handle the fill, then the fill is done via Grid_bcApplyToRegion.
!!
!!  If the given multifab data is defined with respect to a face-centered index
!!  space, then face data in the multifab that coincides with the domain
!!  boundaries will be passed to Grid_bcApplyToRegionSpecialized and 
!!  Grid_bcApplyToRegion.  These routines are allowed to overwrite the data
!!  for those faces that lie on the domain boundary.
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

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine gr_fillPhysicalBC(pmf, scomp, ncomp, time, pgeom) bind(c)
    use iso_c_binding
    
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom
    use amrex_box_module,       ONLY : amrex_box, &
                                       amrex_intersection
    use amrex_geometry_module,  ONLY : amrex_geometry, &
                                       amrex_pmask
    use amrex_multifab_module,  ONLY : amrex_multifab, &
                                       amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy

    use Driver_interface,       ONLY : Driver_abortFlash
    use Grid_data,              ONLY : gr_maxRefine, &
                                       gr_domainBC
    use Grid_interface,         ONLY : Grid_bcApplyToRegion, &
                                       Grid_bcApplyToRegionSpecialized
    use gr_amrexInterface,      ONLY : gr_copyFabInteriorToRegion, &
                                       gr_copyGuardcellRegionToFab

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
   
    type(amrex_box) :: goodData
    type(amrex_box) :: nextGoodData
    type(amrex_box) :: domain
    type(amrex_box) :: guardcells
    type(amrex_box) :: bcData

    integer :: n_cells_domain
    integer :: n_cells_level

    integer :: j
    real    :: delta, delta_j
    integer :: level, face, dir
    integer :: finest_level
    integer :: axis, axis2, axis3

    logical :: mask(ncomp)
    logical :: applied

    integer :: offset

    integer                       :: endPts(LOW:HIGH, 1:MDIM)
    integer                       :: regionType(1:MDIM)
    integer                       :: regionSize(4)
    real(wp), pointer, contiguous :: regionData(:, :, :, :)
    
    real(wp), pointer, contiguous :: solnData(:, :, :, :)

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
    else if (      (.NOT. ntype(IAXIS)) &
             .AND.        ntype(JAXIS) &
             .AND. (.NOT. ntype(KAXIS))) then
        gds = FACEY
    else if (      (.NOT. ntype(IAXIS))  &
             .AND. (.NOT. ntype(JAXIS)) &
             .AND.        ntype(KAXIS)) then
        gds = FACEZ
    else
        call Driver_abortFlash("[gr_fillPhysicalBC] " // &
                               "Given mfab must be cell- or face-centered")
    end if
   
    ! Obtain domain as a box defined w.r.t. the index space of
    ! the given region.  Note that this box contains the domain
    ! boundaries as well.
    !
    ! 0-based, global indices, index space of mfab
    domain = geom%domain
    call domain%convert(ntype)

    call amrex_mfiter_build(mfi, mfab, tiling=.false.)
    do while(mfi%next())
       ! 0-based, global indices, index space of mfab
       solnData => mfab%dataPtr(mfi)

       ! 0-based, global indices, index space of mfab
       box = mfi%fabbox()

       ! Create a box that keeps track across all (face, axis) iterations
       ! that region of the current box that contains correct data.  See above
       ! documentation for what data is assumed correct at the start.
       !
       ! 0-based, global indices, index space of mfab
       goodData = domain
       do axis = 1, NDIM
         if (amrex_pmask(axis)) then
           goodData%lo(axis) = goodData%lo(axis) - NGUARD
           goodData%hi(axis) = goodData%hi(axis) + NGUARD
         end if
       end do
       goodData = amrex_intersection(box, goodData)

       ! A version of goodData that is grown before each step to include those
       ! guardcells in box that will be filled during that step
       nextGoodData = goodData

#ifdef DEBUG_GRID
       write(*,*) "------------------------------------------------------------"
       write(*,*) "------------------------------------------------------------"
       if (gds == CENTER) then
         write(*,*) "CENTER"
       else if (gds == FACEX) then
         write(*,*) "FACEX"
       else if (gds == FACEY) then
         write(*,*) "FACEY"
       else if (gds == FACEZ) then
         write(*,*) "FACEZ"
       end if
       write(*,*) "Domain     Lower = ", domain%lo
       write(*,*) "Domain     Upper = ", domain%hi
       write(*,*) "Box        Lower = ", box%lo
       write(*,*) "Box        Upper = ", box%hi
#endif

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

          offset = 1
          if (     ((gds == FACEX) .AND. (axis == IAXIS)) &
              .OR. ((gds == FACEY) .AND. (axis == JAXIS)) &
              .OR. ((gds == FACEZ) .AND. (axis == KAXIS))) then
              ! Allow the BC handling routine to overwrite face-centered data
              ! that lies on the domain boundary
              offset = 0
          end if

          do face = LOW, HIGH
#ifdef DEBUG_GRID
             write(*,*) "------------------------------------------------------------"
             write(*,*) "Face = ", face, "Axis = ", axis
#endif
             if (     ((face == LOW)  .AND. (box%lo(axis) >= domain%lo(axis))) &
                 .OR. ((face == HIGH) .AND. (box%hi(axis) <= domain%hi(axis))) ) then
#ifdef DEBUG_GRID
                write(*,*) "Doesn't span boundary"
#endif
                CYCLE
             end if

             ! We have configured AMReX to handle periodic BC
             if (gr_domainBC(face, axis) == PERIODIC) then
#ifdef DEBUG_GRID
                write(*,*) "Periodic - AMReX"
#endif
                CYCLE
             end if

             ! We create two blocks here based on goodData that are constructed
             ! w.r.t. to the current domain boundary
             ! - bcData - box of NGUARD layers of current correct data and 
             !            NGUARD layers of guardcells outside of domain boundary
             !            to be filled
             ! - guardcells - the region of guardcells in bcData that are also
             !                in the box that we are filling for AMReX
             !
             ! Both are 0-based, global indices, index space of mfab
             bcData = goodData
             guardcells= goodData
             if (face == LOW) then
                nextGoodData%lo(axis) = box%lo(axis)
                bcData%lo(axis) = goodData%lo(axis) - NGUARD
                bcData%hi(axis) = goodData%lo(axis) + (NGUARD - offset)
                guardcells%lo(axis) = box%lo(axis)
                guardcells%hi(axis) = goodData%lo(axis) - offset
             else
                nextGoodData%hi(axis) = box%hi(axis)
                bcData%lo(axis) = goodData%hi(axis) - (NGUARD - offset)
                bcData%hi(axis) = goodData%hi(axis) + NGUARD
                guardcells%hi(axis) = box%hi(axis)
                guardcells%lo(axis) = goodData%hi(axis) + offset
             end if

#ifdef DEBUG_GRID
             write(*,*) "goodData   Lower = ", goodData%lo
             write(*,*) "goodData   Upper = ", goodData%hi
             write(*,*) "bcData     Lower = ", bcData%lo
             write(*,*) "bcData     Upper = ", bcData%hi
             write(*,*) "guardcells Lower = ", guardcells%lo
             write(*,*) "guardcells Upper = ", guardcells%hi
             write(*,*) "Number GCs = ", guardcells%numpts()
#endif

             ! Convert bcData to FLASH format
             ! 1-based, global indices, index space of mfab
             endPts(:, :) = 1
             endPts(LOW,  1:NDIM) = bcData%lo(1:NDIM) + 1
             endPts(HIGH, 1:NDIM) = bcData%hi(1:NDIM) + 1

             ! Create buffer to hold data with indices permuted for use with
             ! Grid_bcApplyToRegion
             regionSize(BC_DIR)     = endPts(HIGH, axis)  - endPts(LOW, axis)  + 1
             regionSize(SECOND_DIR) = endPts(HIGH, axis2) - endPts(LOW, axis2) + 1
             regionSize(THIRD_DIR)  = endPts(HIGH, axis3) - endPts(LOW, axis3) + 1
             regionSize(STRUCTSIZE) = ncomp

             ! 1-based, local indices, index space of mfab
             allocate(regionData(regionSize(BC_DIR), &
                                 regionSize(SECOND_DIR), &
                                 regionSize(THIRD_DIR), &
                                 regionSize(STRUCTSIZE)) )

             regionData(:, :, :, :) = 0.0
             call gr_copyFabInteriorToRegion(solnData, gds, face, axis, &
                                             goodData, scomp, ncomp, regionData)

             ! As regionData only contains those physical quantities that AMReX
             ! asks for, no need for masking
             mask = .TRUE.

             ! Let simulation do BC fill if so desired
             applied = .FALSE.
             call Grid_bcApplyToRegionSpecialized(gr_domainBC(face, axis), &
                                                  gds, level, NGUARD, &
                                                  axis, face, &
                                                  regionData, regionSize, &
                                                  mask, applied, &
                                                  axis2, axis3, endPts, 0)

             if (.NOT. applied) then
                ! Have FLASH fill GC in special data buffer
                call Grid_bcApplyToRegion(gr_domainBC(face, axis), &
                                          gds, level, NGUARD, &
                                          axis, face, &
                                          regionData, regionSize, &
                                          mask, applied, &
                                          axis2, axis3, endPts, 0)
             end if

             call gr_copyGuardcellRegionToFab(regionData, gds, face, axis, &
                                              guardcells, scomp, ncomp, solnData)

             deallocate(regionData)

             if (.NOT. applied) then
                call Driver_abortFlash("[gr_fillPhysicalBC] BC not applied")
             end if

             goodData = nextGoodData
          end do
       end do

       nullify(solnData)
    end do

    call amrex_mfiter_destroy(mfi)
end subroutine gr_fillPhysicalBC

