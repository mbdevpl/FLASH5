!!****if* source/Grid/localAPI/gr_forceDerefInRadius
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
!!  This implementation is the stub that does nothing.
!!  A real implementation can be found under GridMain/AMR/paramesh .
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

  implicit none

! Arguments

  real, intent(IN)      :: ic, jc, kc, radius
  integer, intent(IN)   :: lref

end subroutine gr_forceDerefInRadius
