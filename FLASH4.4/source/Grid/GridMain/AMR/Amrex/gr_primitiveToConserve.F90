!!****if* source/Grid/GridMain/AMR/Amrex/gr_primitiveToConserve
!!
!! NAME
!!  gr_primitiveToConserve
!!
!! SYNOPSIS
!!  gr_primitiveToConserve(integer(IN) :: lo(MDIM),
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
!!  Given a block of data, convert all cell-centered quantities normally 
!!  represented in mass-specific form (e.g. velocity) into their corresponding
!!  conservative form (e.g. momentum density).  Conversion is achieved by
!!  multiplying primitive-form data by density.  Note that for proper
!!  functioning, DENS_VAR must *not* be marked as PER_MASS.  If density is not
!!  a physical quantity for a simulation, then this conversion is not done.
!!
!!  Cell-centered quantities are considered to be in mass-specific form if they
!!  are explicitly marked as PER_MASS in the Config file.  Additionally,
!!  abundances and mass scalars are considered to be mass-specific.
!!
!!  Note that the inverse conversion requires that density be non-zero on 
!!  all cells.  This routine assumes that all physics units are ensuring that
!!  density is never zero in any cell.
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
!!  gr_conserveToPrimitive
!!
!! BUGS
!!  This routine does not set the variable attributes to indicate that the
!!  variables are now conserved.
!!
!!  This routine accesses the global variable storage array unk directly.  It
!!  won't work for data stored in the workspace array WORK. This method is not
!!  intended for use with a Uniform Grid as its functionality is not needed in
!!  this case.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_primitiveToConserve(lo, hi, &
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

  if (.NOT. ANY(gr_vartypes == VARTYPE_PER_MASS))    RETURN

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Driver_abortFlash('[gr_primitiveToConserve] density is PER_MASS')
  end if

  ! Zero density is non-physical for FLASH simulations and is 
  ! incompatible with the conservative-to-primitive form conversion
  ! done post-interpolation.
  !
  ! Insist that physics units prepare data as they see fit so that
  ! interpolation proceeds with clean data.
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,DENS_VAR) == 0.0) then
          call Driver_abortFlash("[gr_primitiveToConserve] Density is zero")
        end if
      end do
    end do
  end do

  do var = scomp, (scomp + ncomp - 1)
    if (gr_vartypes(var) == VARTYPE_PER_MASS) then
      do     k = lo(KAXIS), hi(KAXIS) 
        do   j = lo(JAXIS), hi(JAXIS) 
          do i = lo(IAXIS), hi(IAXIS)
            d(i,j,k,var) = d(i,j,k,DENS_VAR) * d(i,j,k,var)
          end do
        end do
      end do
    end if
  end do
#endif

end subroutine gr_primitiveToConserve

