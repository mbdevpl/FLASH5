!!****if* source/Grid/GridMain/paramesh/gr_forceDerefInRadius
!!
!! NAME
!!  gr_forceDerefInRadius
!!
!!
!! SYNOPSIS
!!  call gr_forceDerefInRadius(real(in) :: ic,
!!                             real(in) :: jc,
!!                             real(in) :: kc,
!!                             real(in) :: radius,
!!                             integer(in) :: lref)
!!
!! PURPOSE
!!  Mark blocks what lie fully within a certain radius from some point so that
!!  they will be derefined. Either blocks are brought
!!  down to a specific level of refinement or each qualifying block is derefined.
!!
!! ARGUMENTS
!!  ic -   Center of the interval/circle/sphere : IAXIS
!!  jc -                                          JAXIS
!!  kc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius -       Radius of the region
!!  lref  -        If > 0, try to bring all qualifying blocks to this level of refinement;
!!                         if they are already at (exactly) this level, cancel a refine flag;
!!                         if they are at a coarser level, do nothing (allowing refinement if
!!                          refine flag has been set by other parts of the refinement logic).
!!                 If <= 0, derefine qualifying blocks once.
!!
!! NOTES
!!
!!  This routine has been tried only in 1D spherical geometry.
!!
!!  "Forcing" is not absolute - derefinment can still be overridden by PARAMESH
!!  routines if that is required by the criteria for a well-formed grid. In
!!  particular, high refinement levels of blocks that are neighbors to
!!  qualifying blocks can effectively cancel derefine flags and result in
!!  refinement above the desired level.
!!***

subroutine gr_forceDerefInRadius(ic, jc, kc, radius, lref)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  real, intent(IN)      :: ic, jc, kc, radius
  integer, intent(IN)   :: lref

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: bxl, bxr, byl, byr, bzl, bzr
  real                  :: dist2, xdist2, ydist2, zdist2
  integer               :: b


  if((gr_geometry == CARTESIAN).or.(gr_geometry == CYLINDRICAL)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)

           bxl = blockCenter(1) - blockSize(1) - ic
           bxr = blockCenter(1) + blockSize(1) - ic
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
              bzr = 0
           endif

! Find maximum distance from (ic,jc,kc) for each dimension.
! Nonexistent dimensions have had all distances set to zero, above, so
! they are ignored.

           xdist2 = max( bxl**2, bxr**2 )
           ydist2 = max( byl**2, byr**2 )
           zdist2 = max( bzl**2, bzr**2 )

! Now compute the maximum distance to (ic,jc,kc) and compare it to the
! specified radius.  If it is larger than this radius, then the block lies
! at least partially outside the interval/circle/sphere and thus does not
! qualify for derefinement.

           dist2 = xdist2 + ydist2 + zdist2    ! Currently assumes Cartesian
           ! or 2D axisymmetric (r-z)
           ! or 1D spherical (r)
           if (dist2 < radius**2) then ! qualifying block

              if (lrefine(b) > lref ) then
                 derefine(b) = .true.
                 refine(b)   = .false.
              else if (lrefine(b) == lref) then
                 refine(b)   = .false.
              else if (lref <= 0) then
                 derefine(b) = .true.
                 refine(b)   = .false.
              endif

           endif

           ! End of leaf-node block loop
        endif
     end do
  elseif((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL)) then

     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)

           bxl = blockCenter(1) - blockSize(1) - ic
           bxr = blockCenter(1) + blockSize(1) - ic

           dist2 = max( bxl, bxr )

           if (dist2 < radius) then ! qualifying block

              if (lrefine(b) > lref ) then
                 derefine(b) = .true.
                 refine(b)   = .false.
              else if (lrefine(b) == lref) then
                 refine(b)   = .false.
              else if (lref <= 0) then
                 derefine(b) = .true.
                 refine(b)   = .false.
              endif

           endif

        endif
     end do
  else
     call Driver_abortFlash("gr_forceDerefInRadius: geometry spec is wrong")
     !-------------------------------------------------------------------------------
  end if
  return
end subroutine gr_forceDerefInRadius
