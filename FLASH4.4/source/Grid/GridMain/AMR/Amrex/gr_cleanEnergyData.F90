!!****if* source/Grid/GridMain/AMR/Amrex/gr_cleanDensityData
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
!!  This routine performs conservative-to-primitive form conversions.  Refer to
!!  the documentation for gr_conserveToPrimitive for more information.
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
!!  gr_preinterpolationWork
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_cleanEnergyData(smallE, &
                              lo, hi, &
                              d, dlo, dhi, nd)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data,        ONLY : gr_sanitizeDataMode

  implicit none
  
  real,    intent(in)    :: smallE
  integer, intent(in)    :: lo(MDIM), hi(MDIM)
  integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
  integer, intent(in)    :: nd
  real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                              dlo(JAXIS):dhi(JAXIS), &
                              dlo(KAXIS):dhi(KAXIS), &
                              nd)

  integer :: i, j, k

900 format("         Value at",I3,",",I3,",",I3," set to smalle")
#ifdef ENER_VAR
  ! DEV: TODO Determine how to implement all modes and all levels of 
  !           verbosity
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,ENER_VAR) < smallE*0.999999999) then
           ! The factor 0.999999999 is used to avoid triggering aborts
           ! when some physics code, in particular in the Eos unit,
           ! has put a floor of smallE under a very low internal
           ! energy and this then becomes even slightly smaller as a
           ! result of rounding. If it is intended that situations
           ! like that produce messages or aborts, it should be the
           ! responsibility of the physics unit that applied the
           ! smallE floor to take care that that happens.
          if      (gr_sanitizeDataMode == 3) then
            write(*,*) "WARNING: [gr_cleanEnergyData]"
            write(*,*) "         Total energy data less than smalle:", d(i,j,k,ENER_VAR)
            write(*,900)         i, j, k
            d(i,j,k,ENER_VAR) = max(d(i,j,k,ENER_VAR), smallE)
          else if (gr_sanitizeDataMode == 4) then
            call Driver_abortFlash("[gr_cleanEnergyData] ENER data less than smalle")
          end if
        end if
      end do
    end do
  end do
#endif
#ifdef EINT_VAR
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,EINT_VAR) < smallE*0.999999999) then
           ! For the factor 0.999999999 see above.
          if      (gr_sanitizeDataMode == 3) then
            write(*,*) "WARNING: [gr_cleanEnergyData]"
            write(*,*) "         Internal energy data less than smalle:", d(i,j,k,EINT_VAR)
            write(*,900)         i, j, k
            d(i,j,k,EINT_VAR) = max(d(i,j,k,EINT_VAR), smallE)
          else if (gr_sanitizeDataMode == 4) then
            call Driver_abortFlash("[gr_cleanEnergyData] EINT data less than smalle")
          end if
        end if
      end do
    end do
  end do
#endif

end subroutine gr_cleanEnergyData

