!!****if* source/Grid/GridMain/AMR/Amrex/gr_conserveToPrimitive
!!
!! NAME
!!  gr_conserveToPrimitive
!!
!! SYNOPSIS
!!  gr_conserveToPrimitive(block_metadata_t(in) :: block,
!!                         logical(in)          :: allCells)
!!
!! DESCRIPTION
!!  Given a block, convert cell-centered variables that are normally
!!  represented in primitive/mass-specific form (e.g., velocity) from the 
!!  corresponding conservative form (i.e., momentum density) back to the normal
!!  primitive form.  Conversion is achieved by dividing the conservative form by
!!  the density.  Note that for proper functioning, DENS_VAR must *not* be 
!!  marked as PER_MASS.
!!
!!  Cell-centered quantities are considered to be in mass-specific form if they
!!  are explicitly marked as PER_MASS in the Config file.  Additionally,
!!  abundances and mass scalars are considered to be mass-specific.
!!
!!  This conversion is made only if gr_convertToConsvdForMeshCalls is .TRUE.
!!  Otherwise, nothing is done.
!!
!!  A side effect of this routine is that the lower cutoff smalle is applied
!!  to EINT_VAR and ENER_VAR after conversion.  In addition, the lower cutoff
!!  smlrho is applied to DENS_VAR before conversion as for FLASH simulations,
!!  formation of regions of vacuum is considered to be unphysical.  These two
!!  cutoff values are runtime parameters.
!!
!! ARGUMENTS
!!   block - the metadata representation of block whose data shall be
!!           transformed
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
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

subroutine gr_conserveToPrimitive(block, allCells)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,        ONLY : gr_smallrho, &
                               gr_smalle, & 
                               gr_vartypes, &
                               gr_convertToConsvdForMeshCalls
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN) :: block
  logical,                intent(IN) :: allCells
  
  real, pointer :: solnData(:,:,:,:) => null()

  real    :: dens_old_inv
  logical :: needToConvert
  integer :: i, j, k, var
  integer :: ilo, ihi
  integer :: jlo, jhi
  integer :: klo, khi

  if (.NOT. gr_convertToConsvdForMeshCalls)          RETURN

  call Grid_getBlkPtr(block, solnData, CENTER)

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Driver_abortFlash('[gr_conserveToPrimitive] density is PER_MASS')
  end if
  
  if (allCells) then
    ilo = block%limitsGC(LOW,  IAXIS)
    ihi = block%limitsGC(HIGH, IAXIS)
    jlo = block%limitsGC(LOW,  JAXIS)
    jhi = block%limitsGC(HIGH, JAXIS)
    klo = block%limitsGC(LOW,  KAXIS)
    khi = block%limitsGC(HIGH, KAXIS)
  else
    ilo = block%limits(LOW,  IAXIS)
    ihi = block%limits(HIGH, IAXIS)
    jlo = block%limits(LOW,  JAXIS)
    jhi = block%limits(HIGH, JAXIS)
    klo = block%limits(LOW,  KAXIS)
    khi = block%limits(HIGH, KAXIS)
  end if

  needToConvert = ANY(gr_vartypes == VARTYPE_PER_MASS)

  do     k = klo, khi
    do   j = jlo, jhi
      do i = ilo, ihi
        ! cutoff -- in case the interpolants are not monotonic
        solnData(i,j,k,DENS_VAR) = max(solnData(i,j,k,DENS_VAR), gr_smallrho)

        if (needToConvert) then
          dens_old_inv = 1.0 / solnData(i,j,k,DENS_VAR)

          do var = UNK_VARS_BEGIN, UNK_VARS_END
            if (gr_vartypes(var) == VARTYPE_PER_MASS) then
              solnData(i,j,k,var) = dens_old_inv * solnData(i,j,k,var)
            end if
          end do
        end if

      end do
    end do
  end do
#endif

#ifdef ENER_VAR               
  solnData(:,:,:,ENER_VAR) = max(solnData(:,:,:,ENER_VAR), gr_smalle)
#endif
#ifdef EINT_VAR
  solnData(:,:,:,EINT_VAR) = max(solnData(:,:,:,EINT_VAR), gr_smalle)
#endif

  call Grid_releaseBlkPtr(block, solnData, CENTER)

end subroutine gr_conserveToPrimitive
 
