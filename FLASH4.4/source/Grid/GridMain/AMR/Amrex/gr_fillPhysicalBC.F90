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
!!
!!
!! ARGUMENTS 
!!  
!!
!!
!! NOTES
!!
!!
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
    use Grid_interface,         ONLY : Grid_bcApplyToRegion
    use gr_amrexInterface,      ONLY : gr_getPatchBoundaryEndpoints, &
                                       gr_transformBcRegion, &
                                       gr_untransformBcRegion
    use gr_physicalMultifabs,   ONLY : unk
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

    real    :: delta, delta_j
    integer :: i, j, k, var
    integer :: level, face, dir
    integer :: finest_level
    integer :: axis, axis2, axis3

    integer :: intersection(LOW:HIGH, 1:MDIM)
    logical :: mask(NUNK_VARS)
    logical :: applied

    integer                       :: endPts(LOW:HIGH, 1:MDIM)
    integer                       :: regionType(1:MDIM)
    integer                       :: regionSize(4)
    real(wp), pointer, contiguous :: regionData(:, :, :, :)
    
    real(wp), pointer, contiguous :: solnData(:, :, :, :)
    integer                       :: limitsGC(LOW:HIGH, 1:MDIM)
    
    if (amrex_is_all_periodic())    RETURN

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
       ! 0-based, cell-centered, and global
       box = mfi%fabbox()
       solnData => mfab%dataPtr(mfi)

       ! 0-based, cell-centered, and global indices
       limitsGC(:, :) = 0
       limitsGC(LOW,  1:NDIM) = box%lo(1:NDIM)
       limitsGC(HIGH, 1:NDIM) = box%hi(1:NDIM)

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
             ! Get end points that define region that needs BC filling
             ! endPts is 1-based, cell-centered, and global
             call gr_getPatchBoundaryEndpoints(face, axis, limitsGC, &
                                               amrex_geom(level-1)%dx, endpts)

             ! Skip if face not on boundary
             if (endpts(LOW, axis) == endpts(HIGH, axis))   CYCLE

             intersection(:, :) = 0
             do dir = 1, NDIM
                intersection(LOW,  dir) = MAX(limitsGC(LOW,  dir), &
                                              endPts(LOW,  dir)-1)
                intersection(HIGH, dir) = MIN(limitsGC(HIGH, dir), &
                                              endPts(HIGH, dir)-1)
             end do

             ! Create buffer to hold data with indices permuted for use with
             ! Grid_bcApplyToRegion
             regionSize(BC_DIR)     = endpts(HIGH, axis)  - endpts(LOW, axis)  + 1
             regionSize(SECOND_DIR) = endpts(HIGH, axis2) - endpts(LOW, axis2) + 1
             regionSize(THIRD_DIR)  = endpts(HIGH, axis3) - endpts(LOW, axis3) + 1
             regionSize(STRUCTSIZE) = NUNK_VARS

             ! 1-based, cell-centered, and LOCAL
             allocate(regionData(regionSize(BC_DIR), &
                                 regionSize(SECOND_DIR), &
                                 regionSize(THIRD_DIR), &
                                 regionSize(STRUCTSIZE)) )

             ! Populate regionData with given interior data
             call gr_transformBcRegion(solnData, axis, intersection, &
                                       regionSize, regionData)

             ! Fill buffer with BC
             mask = .TRUE.
             applied = .FALSE.
             call Grid_bcApplyToRegion(gr_domainBC(face, axis), CENTER, NGUARD, &
                                       axis, face, regionData, regionSize, &
                                       mask, applied, blockDesc, &
                                       axis2, axis3, endpts, 0)

             ! Copy data from buffer to target
             call gr_untransformBcRegion(regionData, axis, intersection, &
                                         regionSize, solnData)

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

