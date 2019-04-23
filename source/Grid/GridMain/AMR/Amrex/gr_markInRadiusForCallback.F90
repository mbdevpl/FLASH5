!!****if* source/Grid/GridMain/AMR/Amrex/gr_markInRadiusForCallback
!!
!! NAME
!!  gr_markInRadius
!!
!!
!! SYNOPSIS
!!  gr_markInRadius(real(in) :: ic,
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

subroutine gr_markInRadiusForCallback(ic, jc, kc, radius, lev, tags, tagval)

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

  real,                    intent(IN) :: ic, jc, kc, radius
  integer,                 intent(IN) :: lev
  type(c_ptr),             intent(IN) :: tags
  character(c_char),       intent(IN) :: tagval

! Local data

  type(amrex_tagboxarray) :: tag
  type(amrex_mfiter)      :: mfi
  type(amrex_box)         :: bx
  type(Grid_tile_t)       :: blockDesc

  character(c_char), contiguous, pointer :: tagData(:,:,:,:)

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: bxl, bxr, byl, byr, bzl, bzr
  real                  :: dist2, xdist2, ydist2, zdist2

  integer :: i, j, k

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
     blockSize(:)  = 0.5 * blockSize(:)

! Find minimum distance from (ic,jc,kc) for each dimension.  For each
! coordinate, if both "left" and "right" distances have the same sign,
! then the smaller magnitude is the minimum.  Otherwise (ic,jc,kc) is
! contained within the interval for that dimension, so the minimum is 0.
! Nonexistent dimensions have had all distances set to zero, so they are
! ignored.

     bxl = blockCenter(1) - blockSize(1) - ic
     bxr = blockCenter(1) + blockSize(1) - ic

     if ((gr_geometry == CARTESIAN).or.(gr_geometry == CYLINDRICAL)) then
        if (NDIM > 1) then
           byl = blockCenter(2) - blockSize(2) - jc
           byr = blockCenter(2) + blockSize(2) - jc
        else
           byl = 0.
           byr = 0.
        endif
        if ((NDIM == 3).and.(gr_geometry==CARTESIAN)) then
           bzl = blockCenter(3) - blockSize(3) - kc
           bzr = blockCenter(3) + blockSize(3) - kc
        else
           bzl = 0.
           bzr = 0.
        endif

! Now compute the minimum distance to (ic,jc,kc) and compare it to the
! specified radius.  If it is less than this radius, then the block contains
! at least part of the interval/circle/sphere and is marked for refinement.

        if (bxl*bxr > 0.) then
           xdist2 = min( bxl**2, bxr**2 )
        else
           xdist2 = 0.
        endif
        if (byl*byr > 0.) then
           ydist2 = min( byl**2, byr**2 )
        else
           ydist2 = 0.
        endif
        if (bzl*bzr > 0.) then
           zdist2 = min( bzl**2, bzr**2 )
        else
           zdist2 = 0.
        endif
        dist2 = xdist2 + ydist2 + zdist2
     elseif ((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL)) then
        if (bxl*bxr > 0.) then
           dist2 = min( bxl**2, bxr**2 )
        else
           dist2 = 0.
        endif
     else
        call Driver_abortFlash("MarkRefine: geometry spec is wrong")
     endif

     tagData => tag%dataptr(mfi)

     associate (lo     => blockDesc%limits(LOW,  :), &
                hi     => blockDesc%limits(HIGH, :), &
                lo_tag => lbound(tagData), &
                hi_tag => ubound(tagData))

        if (dist2 <= radius**2) then
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
end subroutine gr_markInRadiusForCallback
