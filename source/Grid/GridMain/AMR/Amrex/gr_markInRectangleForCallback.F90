!!****if* source/Grid/GridMain/AMR/Amrex/gr_markInRectangleForCallback
!!
!! NAME
!!  gr_markInRectangle
!!
!!
!! SYNOPSIS
!!  gr_markInRectangle(real(in) :: ic,
!!                  real(in) :: jc,
!!                  real(in) :: kc,
!!                  real(in) :: radius,
!!                  integer(in) :: lev,
!!                  c_ptr(in) :: tags,
!!                  c_char(in) :: tagval)
!!
!! PURPOSE
!!  Refine all blocks containing points within a circular/spherical region of
!!  given radius about a given point (xc,yc,zc).  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.
!!
!! ARGUMENTS
!!  ic -   Center of the interval/circle/sphere : IAXIS
!!  jc -                                          JAXIS
!!  kc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius -       Radius of the region
!!  lev - the 0-based level index
!!  tags - C-pointer to an AMReX tagbox array.  The elements of this are tag
!!         boxes.  The cells of these tagboxes are set to communicate a need
!!         to refine the associated block.
!!  tagval - for full, rich AMReX tagging, this values should be assigned to
!!           each cell that has insufficient resolution.
!!
!! NOTES
!!
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!
!!***

subroutine gr_markInRectangleForCallback(ilb, irb, jlb, jrb, klb, krb, contained, &
                                         lev, tags, tagval)

!-------------------------------------------------------------------------------
  use iso_c_binding
  use amrex_fort_module,      ONLY : wp => amrex_real
  use amrex_box_module,       ONLY : amrex_box
  use amrex_tagbox_module,    ONLY : amrex_tagboxarray
  use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                     amrex_mfiter_build, &
                                     amrex_mfiter_destroy

  use Driver_interface,       ONLY : Driver_abortFlash
  use Grid_data,              ONLY : gr_geometry
  use Grid_interface,         ONLY : Grid_getBlkCenterCoords
  use gr_physicalMultifabs,   ONLY : unk
  use Grid_tile,              ONLY : Grid_tile_t
#include "constants.h"
#include "Flash.h"

  implicit none

! Arguments

  real,                    intent(IN) :: ilb, irb, jlb, jrb, klb, krb
  integer,                 intent(IN) :: lev, contained
  type(c_ptr),             intent(IN) :: tags
  character(c_char),       intent(IN) :: tagval

! Local data

  type(amrex_tagboxarray) :: tag
  type(amrex_mfiter)      :: mfi
  type(amrex_box)         :: bx
  type(Grid_tile_t)       :: blockDesc

  character(c_char), contiguous, pointer :: tagData(:,:,:,:)

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: xl, xr, yl, yr, zl, zr
  integer               :: b
  logical               :: x_in_rect, y_in_rect, z_in_rect

  integer :: i, j, k

#ifdef DEBUG
  if((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL))&
       call Driver_abortFlash("markRefineInRectangle : wrong geometry")
  if((gr_geometry==CYLINDRICAL).and.(NDIM==3))&
       call Driver_abortFlash("markRefineInRectangle : not valid in 3d for cylindrical")
#endif
  tag = tags

  call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)
  do while(mfi%next())
     bx = mfi%fabbox()

     blockDesc%level = lev + 1
     blockDesc%grid_index = mfi%grid_index()
     blockDesc%limits(LOW,  :) = 1
     blockDesc%limits(HIGH, :) = 1
     blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1 + NGUARD
     blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1 - NGUARD
     blockDesc%blkLimitsGC(LOW,  :) = 1
     blockDesc%blkLimitsGC(HIGH, :) = 1
     blockDesc%blkLimitsGC(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
     blockDesc%blkLimitsGC(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
     blockDesc%grownLimits(:, :) = blockDesc%blkLimitsGC(:, :)

     call Grid_getBlkCenterCoords(blockDesc, blockCenter)
     call blockDesc % physicalSize(blockSize)

        xl = blockCenter(1) - blockSize(1)
        xr = blockCenter(1) + blockSize(1)
        if (NDIM > 1) then
           yl = blockCenter(2) - blockSize(2)
           yr = blockCenter(2) + blockSize(2)
        endif
        if (NDIM == 3) then
           zl = blockCenter(3) - blockSize(3)
           zr = blockCenter(3) + blockSize(3)
        endif

        ! For each dimension, determine whether the block overlaps the specified
        ! rectangle.  Nonexistent dimensions are ignored.  This method assumes
        ! Cartesian coordinates (or the cross-section of a rectangular torus in
        ! 2D axisymmetric coordinates, or an annulus in 1D spherical coordinates).

        if (contained /= 0) then       ! only refine if completely contained

           x_in_rect =   ((xl >= ilb) .and. (xr <= irb))

           if (NDIM >= 2) then
              y_in_rect = ((yl >= jlb) .and. (yr <= jrb))
           else
              y_in_rect = .true.
           endif

           if (NDIM == 3) then
              z_in_rect = ((zl >= klb) .and. (zr <= krb))
           else
              z_in_rect = .true.
           endif

        else                           ! refine if any overlap with rectangle

           x_in_rect =   .not. ((xr <= ilb) .or. (xl >= irb))

           if (NDIM >= 2) then
              y_in_rect = .not. ((yr <= jlb) .or. (yl >= jrb))
           else
              y_in_rect = .true.
           endif

           if (NDIM == 3) then
              z_in_rect = .not. ((zr <= klb) .or. (zl >= krb))
           else
              z_in_rect = .true.
           endif

        endif

     ! Refine the block if all of the dimensions overlap/are contained.

     tagData => tag%dataptr(mfi)

     associate (lo     => blockDesc%limits(LOW,  :), &
                hi     => blockDesc%limits(HIGH, :), &
                lo_tag => lbound(tagData), &
                hi_tag => ubound(tagData))

        if (x_in_rect .and. y_in_rect .and. z_in_rect) then
           i = INT(0.5d0 * DBLE(lo_tag(IAXIS) + hi_tag(IAXIS)))
           j = INT(0.5d0 * DBLE(lo_tag(JAXIS) + hi_tag(JAXIS)))
           k = INT(0.5d0 * DBLE(lo_tag(KAXIS) + hi_tag(KAXIS)))

           ! Fourth index is 1:1
           tagData(i, j, k, 1) = tagval
        endif

     end associate

     nullify(tagData)
  end do
  call amrex_mfiter_destroy(mfi)

  return
end subroutine gr_markInRectangleForCallback
