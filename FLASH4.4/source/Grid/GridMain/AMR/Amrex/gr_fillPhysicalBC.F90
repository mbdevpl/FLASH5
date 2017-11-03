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
    use Grid_data,              ONLY : gr_maxRefine
    use Grid_interface,         ONLY : Grid_getBlkBC
    use gr_bcInterface,         ONLY : gr_bcApplyToOneFace
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
    integer :: j, level, axis
    integer :: finest_level

    integer :: regionType(1:MDIM)
    integer :: faces(LOW:HIGH, 1:MDIM)
 
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

    ! Let FLASH fill GC at physical boundaries for all blocks in level
    ! DEV: FIXME Using mfab instead of unk here lead to iterating over boxes
    ! that were not really blocks.  What is going on?
    call amrex_mfiter_build(mfi, unk(level-1), tiling=.false.)

    do while(mfi%next())
       ! Manually create block descriptor
       box = mfi%tilebox()

       blockDesc%level = level
       blockDesc%grid_index = mfi%grid_index()
       blockDesc%limits(LOW,  :) = 1
       blockDesc%limits(HIGH, :) = 1
       blockDesc%limits(LOW,  1:NDIM) = box%lo(1:NDIM) + 1
       blockDesc%limits(HIGH, 1:NDIM) = box%hi(1:NDIM) + 1
       blockDesc%limitsGC(LOW,  :) = 1
       blockDesc%limitsGC(HIGH, :) = 1
       blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
       blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

       call Grid_getBlkBC(blockDesc, faces)

!       write(*,*) "[gr_fillPhysicalBC] Block ", blockDesc%grid_index
!       write(*,*) "[gr_fillPhysicalBC] Level = ", blockDesc%level
!       write(*,*) "[gr_fillPhysicalBC] Limits = ", blockDesc%limits(LOW,  :)
!       write(*,*) "[gr_fillPhysicalBC] Limits = ", blockDesc%limits(HIGH, :)
!       write(*,*) "[gr_fillPhysicalBC] Lower Faces = ", faces(LOW,  :)
!       write(*,*) "[gr_fillPhysicalBC] Upper Faces = ", faces(HIGH, :)

       ! DEV: TODO I am passing CENTER here as we are just dealing with unk
       ! at the moment.  Otherwise, we will need to infer the FACE[XYZ] from
       ! the index type of the multifab.  Check if we can even get this
       ! through the Fortran interface yet.

       ! DEV: TODO Handle corner block case where no face is on the physical
       ! but the corner GC region is on the physical domain.

        do axis = 1, NDIM
            if (faces(LOW, axis) /= NOT_BOUNDARY) then
                regionType(:) = WHOLE_VECTOR
                regionType(axis) = LEFT_EDGE
                call gr_bcApplyToOneFace(axis, faces(LOW, axis), &
                                         CENTER, NUNK_VARS, &
                                         regionType, blockDesc, 0)
            end if
            
            if (faces(HIGH, axis) /= NOT_BOUNDARY) then
                regionType(:) = WHOLE_VECTOR
                regionType(axis) = RIGHT_EDGE
                call gr_bcApplyToOneFace(axis, faces(HIGH, axis), &
                                         CENTER, NUNK_VARS, &
                                         regionType, blockDesc, 0)
            end if
        end do
    end do

    call amrex_mfiter_destroy(mfi)

end subroutine gr_fillPhysicalBC

