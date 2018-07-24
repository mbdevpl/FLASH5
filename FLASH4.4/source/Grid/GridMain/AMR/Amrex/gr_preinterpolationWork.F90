!!****if* source/Grid/GridMain/AMR/Amrex/gr_preinterpolationWork
!!
!! NAME
!!  gr_preinterpolationWork
!!
!! SYNOPSIS
!!  gr_preinterpolationWork(integer(IN)       :: lo(MDIM),
!!                          integer(IN)       :: hi(MDIM),
!!                          amrex_real(INOUT) :: d(:,:,:,:),
!!                          integer(IN)       :: dlo(MDIM),
!!                          integer(IN)       :: dhi(MDIM),
!!                          integer(IN)       :: nd,
!!                          integer(IN)       :: scomp,
!!                          integer(IN)       :: ncomp)
!!
!! DESCRIPTION
!!  This is a callback routine that is passed to the AMReX fill patch
!!  routines that may carry out interpolation as part of their data 
!!  movement actions.
!!
!!  AMReX calls this routine just before carrying out interpolation so that
!!  the calling application may perform any necessary work on the data
!!  that is used to perform the interpolation.  This allows for altering
!!  the data only when necessary, rather than altering *all* data before
!!  calling a routine that *might* perform interpolation.
!!
!!  For the GC fill, this routine allows for transforming primitive form data
!!  to conservative form on the interior cells of a coarse block at a 
!!  fine/coarse boundary.
!!
!!  When creating a new leaf block, this routine allows for transforming 
!!  primitive form data to conservative form on the interior and guardcells
!!  of the parent block. 
!!
!!  This routine and its post-interpolation partner together perform
!!  primitive-to-conservative and conservative-to-primitive form conversions.
!!  To accomplish the latter, it is required that the density data be
!!  non-zero.  Therefore, this routine requires that all physics units that 
!!  initiate interpolation must first ensure that all density data is non-zero.
!!
!!  This routine should never be called directly within FLASH.
!!
!! ARGUMENTS
!!  lo - the index of the lower-left corner of the region of the given
!!       box that requires conversion from primitive-to-conservative.
!!  hi - the index of the upper-right corner of the region of the given
!!       box that requires cleaning and conversion from
!!       conservative-to-primitive.
!!  dlo - the lower-left index of the given box
!!  dhi - the upper-right index of the given box
!!  nd -  the number of physical quantities in box
!!  scomp - the index of the first physical quantity in the box
!!  ncomp - the number of physical quantities in the box
!!  d - the data array for the box that requires interpolation
!!
!! SEE ALSO
!!  gr_postinterpolationWork
!!  Grid_fillGuardcells
!!  Grid_updateRefinement
!!  gr_makeFineLevelFromCoarseCallback
!!  gr_remakeLevelCallback
!!
!!***

#include "constants.h"

subroutine gr_preinterpolationWork(lo, hi, &
                                   d, dlo, dhi, nd, &
                                   scomp, ncomp) bind(c)
  use iso_c_binding,     ONLY : c_int

  use amrex_fort_module, ONLY : wp => amrex_real

  use Driver_interface,  ONLY : Driver_abortFlash
  use Grid_data,         ONLY : gr_convertToConsvdInMeshInterp
  use gr_amrexInterface, ONLY : gr_primitiveToConserve

  implicit none

  integer(c_int), intent(in)          :: lo(MDIM), hi(MDIM)
  integer(c_int), intent(in)          :: dlo(MDIM), dhi(MDIM)
  integer(c_int), intent(in),   value :: nd
  integer(c_int), intent(in),   value :: scomp
  integer(c_int), intent(in),   value :: ncomp
  real(wp),       intent(inout)       :: d(dlo(IAXIS):dhi(IAXIS), &
                                           dlo(JAXIS):dhi(JAXIS), &
                                           dlo(KAXIS):dhi(KAXIS), &
                                           nd)

  if (gr_convertToConsvdInMeshInterp) then
    call gr_primitiveToConserve(lo, hi, d, dlo, dhi, nd, scomp, ncomp)
  end if

end subroutine gr_preinterpolationWork 

