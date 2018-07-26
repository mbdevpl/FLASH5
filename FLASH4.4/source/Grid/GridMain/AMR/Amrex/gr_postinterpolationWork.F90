!!****if* source/Grid/GridMain/AMR/Amrex/gr_postinterpolationWork
!!
!! NAME
!!  gr_postinterpolationWork
!!
!! SYNOPSIS
!!  gr_postinterpolationWork(integer(IN)       :: lo(MDIM),
!!                           integer(IN)       :: hi(MDIM),
!!                           amrex_real(INOUT) :: d(dlo(IAXIS):dhi(IAXIS), &
!!                                                  dlo(JAXIS):dhi(JAXIS), &
!!                                                  dlo(KAXIS):dhi(KAXIS), &
!!                                                  nd),
!!                           integer(IN)       :: dlo(MDIM),
!!                           integer(IN)       :: dhi(MDIM),
!!                           integer(IN)       :: nd,
!!                           integer(IN)       :: scomp,
!!                           integer(IN)       :: ncomp)
!!
!! DESCRIPTION
!!  This is a callback routine that is passed to the AMReX fill patch
!!  routines, which may carry out interpolation as part of their data 
!!  movement actions.
!!
!!  AMReX calls this routine just after carrying out interpolation so that
!!  the calling application may perform any necessary work on the data
!!  that was used to perform the interpolation and on the data that was
!!  set via interpolation.
!!
!!  This routine performs conservative-to-primitive form conversion.
!!
!!  Interpolation might lead to non-physical data values.  Before converting, 
!!  all density values less than the runtime parameter smlrho are set to smlrho.
!!  After conversion, the lower cutoff smalle runtime parameter is similarly
!!  applied to EINT_VAR and ENER_VAR if these are quantities in the simulation.
!!
!!  For the GC fill, this routine allows for transforming conservative form data
!!  back to primitive form on the interior cells of a coarse block at a 
!!  fine/coarse boundary.  The interpolated GC data must also be transformed
!!  to primitive form.
!!
!!  When creating a new leaf block, this routine allows for transforming 
!!  conservative form data back to primitive form on the interior and guardcells
!!  of the parent block.  The interpolated interior and GC data of the leaf 
!!  block must also be transformed to primitive form.
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
!!  gr_conserveToPrimitive
!!  gr_cleanDensityData
!!  gr_cleanEnergyData
!!  gr_preinterpolationWork
!!
!!***

#include "constants.h"

subroutine gr_postinterpolationWork(lo, hi, &
                                    d, dlo, dhi, nd, &
                                    scomp, ncomp) bind(c)
  use iso_c_binding,     ONLY : c_int

  use amrex_fort_module, ONLY : wp => amrex_real
  
  use Grid_data,         ONLY : gr_convertToConsvdInMeshInterp, &
                                gr_smallrho, &
                                gr_smalle
  use gr_amrexInterface, ONLY : gr_conserveToPrimitive

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

  ! Zero density values are considered as errors during the
  ! pre-interpolation phase.  Therefore, zero density here, which
  ! is considered to be non-physical, is the result of interpolation.
  ! Hence, we correct before converting to primitive form.
  call gr_cleanDensityData(gr_smallrho, lo, hi, d, dlo, dhi, nd)

  if (gr_convertToConsvdInMeshInterp) then
    call gr_conserveToPrimitive(lo, hi, d, dlo, dhi, nd, scomp, ncomp)
  end if

  ! Zero energy is unphysical for FLASH simulations.  It is assumed that physics
  ! units correct for unacceptable values before initiating interpolation.  
  ! Therefore, incorrect values here must arise from interpolation.  Hence,
  ! we clean these once after they have been potentiall reverted to primitive
  ! form
  call gr_cleanEnergyData(gr_smalle, lo, hi, d, dlo, dhi, nd)

end subroutine gr_postinterpolationWork 

