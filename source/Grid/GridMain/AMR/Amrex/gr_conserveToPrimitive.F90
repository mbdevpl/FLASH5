!!****if* source/Grid/GridMain/AMR/Amrex/gr_conserveToPrimitive
!!
!! NAME
!!  gr_conserveToPrimitive
!!
!! SYNOPSIS
!!  gr_conserveToPrimitive(integer(IN) :: lo(MDIM),
!!                         integer(IN) :: hi(MDIM),
!!                         real(INOUT) :: d(dlo(IAXIS):dhi(IAXIS), &
!!                                          dlo(JAXIS):dhi(JAXIS), &
!!                                          dlo(KAXIS):dhi(KAXIS), &
!!                                          nd),
!!                         integer(IN) :: dlo(MDIM),
!!                         integer(IN) :: dhi(MDIM),
!!                         integer(IN) :: nd,
!!                         integer(IN) :: scomp,
!!                         integer(IN) :: ncomp)
!!
!! DESCRIPTION
!!  Given a block of data, convert cell-centered variables that are normally
!!  represented in primitive/mass-specific form (e.g., velocity) from the 
!!  corresponding conservative form (i.e., momentum density) back to the normal
!!  primitive form.  Conversion is achieved by dividing the conservative form by
!!  the density.  Calling functions must ensure that the density in non-zero
!!  everywhere.  Note that for proper functioning, DENS_VAR must *not* be 
!!  marked as PER_MASS.  If density is note a physical quantity for a
!!  simulation, then this conversion is not done.
!!
!!  Cell-centered quantities are considered to be in mass-specific form if they
!!  are explicitly marked as PER_MASS in the Config file.  Additionally,
!!  abundances and mass scalars are considered to be mass-specific.
!!
!!  This routine does not check if this conversion is enabled by runtime 
!!  parameters.  Rather, the calling routine must know that this conversion
!!  is necessary and desired.
!!
!! ARGUMENTS
!!  lo/hi - the lower and upper corners that define the region of cells
!!          in the given block of data on which the conversion shall be done
!!  dlo/dhi/nd - the lower and upper bounds of the indices of the given data
!!  scomp - the first physical quantity to be potentially converted
!!  ncomp - the number of physical quantities to be potentially converted
!!  d - the data to convert
!!
!! SEE ALSO
!!  gr_primitiveToConserve 
!!
!! BUGS
!!  This routine does not set the variable attributes to indicate that the
!!  variables are no longer conserved.
!!
!!  This routine accesses the global variable storage array unk directly.  It
!!  won't work for data stored in the workspace array WORK. This method is not
!!  intended for use with a Uniform Grid as its functionality is not needed in
!!  this case.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_conserveToPrimitive(lo, hi, &
                                  d, dlo, dhi, nd, &
                                  scomp, ncomp)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data,        ONLY : gr_vartypes

  implicit none

  integer, intent(in)    :: lo(MDIM), hi(MDIM)
  integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
  integer, intent(in)    :: nd
  integer, intent(in)    :: scomp
  integer, intent(in)    :: ncomp
  real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                              dlo(JAXIS):dhi(JAXIS), &
                              dlo(KAXIS):dhi(KAXIS), &
                              nd)

  integer :: i, j, k, var

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Driver_abortFlash('[gr_conserveToPrimitive] density is PER_MASS')
  end if

  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,DENS_VAR) == 0.0) then
          call Driver_abortFlash("[gr_conserveToPrimitive] Density is zero")
        end if
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

end subroutine gr_conserveToPrimitive
 
