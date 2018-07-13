!!****if* source/Grid/GridMain/AMR/Amrex/gr_postinterpolationWork
!!
!! NAME
!!  gr_postinterpolationWork
!!
!! SYNOPSIS
!!  gr_postinterpolationWork(integer(IN)       :: lo(MDIM),
!!                           integer(IN)       :: hi(MDIM),
!!                           amrex_real(INOUT) :: d(:,:,:,:),
!!                           integer(IN)       :: dlo(MDIM),
!!                           integer(IN)       :: dhi(MDIM),
!!                           integer(IN)       :: nd,
!!                           integer(IN)       :: scomp,
!!                           integer(IN)       :: ncomp)
!!
!! DESCRIPTION
!!  This is a callback routine that is passed to the AMReX fill patch
!!  routines that may carry out interpolation as part of their data 
!!  movement actions.
!!
!!  AMReX calls this routine just after carrying out interpolation so that
!!  the calling application may perform any necessary work on the data
!!  that was used to perform the interpolation and on the data that
!!  set via interpolation.
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
!!  This routine performs conservative-to-primitive form conversions,
!!  which requires that the density data be non-zero.  While the
!!  pre-interpolation insists that this be true, it might not be the case
!!  after interpolation.  Therefore, this routine will threshold the density
!!  data so that all data is greater than or equal to the value of the smlrho
!!  runtime parameter.
!!
!!  In addition, if ENER_VAR and EINT_VAR are present in the simulation and if
!!  interpolation is performed on this quantities, then this data will also
!!  be threshold so that all data is greater than or equal to the value of the
!!  smalle runtime parameter.
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
!!  gr_preinterpolationWork
!!  Grid_fillGuardcells
!!  Grid_updateRefinement
!!  gr_makeFineLevelFromCoarseCallback
!!  gr_remakeLevelCallback
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_postinterpolationWork(lo, hi, &
                                    d, dlo, dhi, nd, &
                                    scomp, ncomp) bind(c)
  use iso_c_binding, ONLY : c_int

  use amrex_fort_module, ONLY : wp => amrex_real
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data,        ONLY : gr_convertToConsvdInMeshInterp, &
                               gr_vartypes, &
                               gr_smallrho, &
                               gr_smalle

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

  integer :: i, j, k, var
   
   ! DEV FIXME: This is termporarily commented out so that we can test
   ! conservation as these hooks are phased in
!  if (.NOT. gr_convertToConsvdInMeshInterp)    RETURN

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Driver_abortFlash('[gr_postinterpolationWork] Density is PER_MASS')
  end if

  if ((DENS_VAR < scomp) .OR. (DENS_VAR > scomp+ncomp-1)) then
    call Driver_abortFlash("[gr_postinterpolationWork] Density data not given")
  end if

  ! Zero density values are considered as errors during the
  ! pre-interpolation phase.  Therefore, zero density here, which
  ! is considered to be non-physical, is the result of interpolation.
  ! Hence, we correct before converting to primitive form.
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        d(i,j,k,DENS_VAR) = max(d(i,j,k,DENS_VAR), gr_smallrho)
      end do
    end do
  end do

  do var = scomp, (scomp + ncomp - 1)
    if (gr_vartypes(var) == VARTYPE_PER_MASS) then
      do     k = lo(KAXIS), hi(KAXIS) 
        do   j = lo(JAXIS), hi(JAXIS) 
          do i = lo(IAXIS), hi(IAXIS)
            d(i,j,k,var) = d(i,j,k,var) / d(i,j,k,DENS_VAR)
          end do
        end do
      end do
    end if
  end do
#endif

#ifdef ENER_VAR
  ! Zero energy is unphysical for FLASH simulations.  It is assumed that physics
  ! units correct for unacceptable values before initiating interpolation.  
  ! Therefore, incorrect values here must arise from interpolation.
  if ((ENER_VAR >= scomp) .AND. (ENER_VAR <= scomp+ncomp-1)) then
    do     k = lo(KAXIS), hi(KAXIS) 
      do   j = lo(JAXIS), hi(JAXIS) 
        do i = lo(IAXIS), hi(IAXIS)
          d(i,j,k,ENER_VAR) = max(d(i,j,k,ENER_VAR), gr_smalle)
        end do
      end do
    end do
  end if
#endif

#ifdef EINT_VAR
  ! See comment in ENER_VAR pre-processor block
  if ((EINT_VAR >= scomp) .AND. (EINT_VAR <= scomp+ncomp-1)) then
    do     k = lo(KAXIS), hi(KAXIS) 
      do   j = lo(JAXIS), hi(JAXIS) 
        do i = lo(IAXIS), hi(IAXIS)
          d(i,j,k,EINT_VAR) = max(d(i,j,k,EINT_VAR), gr_smalle)
        end do
      end do
    end do
  end if
#endif

end subroutine gr_postinterpolationWork 

