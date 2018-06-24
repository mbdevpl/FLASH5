!!****if* source/Grid/GridMain/AMR/Amrex/gr_primitiveToConserve
!!
!! NAME
!!  gr_primitiveToConserve
!!
!! SYNOPSIS
!!  gr_primitiveToConserve(block_metadata_t(IN) :: block)
!!
!! DESCRIPTION
!!  Given a block, convert all cell-centered quantities normally represented in
!!  mass-specific form (e.g. velocity) into their corresponding conservative
!!  form (e.g. momentum density).  Conversion is achieved by multiplying
!!  primitive-form data by density.  Note that for proper functioning, DENS_VAR
!!  must *not* be marked as PER_MASS.
!!
!!  Cell-centered quantities are considered to be in mass-specific form if they
!!  are explicitly marked as PER_MASS in the Config file.  Additionally,
!!  abundances and mass scalars are considered to be mass-specific.
!!
!!  This conversion is made only if gr_convertToConsvdForMeshCalls is .TRUE.
!!  Otherwise, nothing is done.
!!
!!  A side effect of this routine is that the lower cutoff smalle is applied
!!  to EINT_VAR and ENER_VAR before conversion.  In addition, the lower cutoff
!!  smlrho is applied to DENS_VAR before conversion as for FLASH simulations,
!!  formation of regions of vacuum is considered to be unphysical.  These two
!!  cutoff values are runtime parameters.
!!
!! ARGUMENTS
!!   block - the metadata representation of block whose data shall be
!!           transformed
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

subroutine gr_primitiveToConserve(block)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,        ONLY : gr_convertToConsvdForMeshCalls, &
                               gr_vartypes, &
                               gr_smalle, &
                               gr_smallrho
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN) :: block

  real, pointer :: solnData(:,:,:,:) => null()
  integer       :: i, j, k, var
  logical       :: fixedDensity

  if (.NOT. ANY(gr_vartypes == VARTYPE_PER_MASS))    RETURN
  if (.NOT. gr_convertToConsvdForMeshCalls)          RETURN

  call Grid_getBlkPtr(block, solnData, CENTER)

  ! Apply energy cutoffs before possible conversion/use
#ifdef ENER_VAR
  solnData(:,:,:,ENER_VAR) = max(solnData(:,:,:,ENER_VAR), gr_smalle)
#endif
#ifdef EINT_VAR
  solnData(:,:,:,EINT_VAR) = max(solnData(:,:,:,EINT_VAR), gr_smalle)
#endif

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Grid_releaseBlkPtr(block, solnData, CENTER)
    call Driver_abortFlash('[gr_primitiveToConserve] density is PER_MASS')
  end if

  ! DEV: TODO Add non-permanent GC here
  associate(ilo => block%limits(LOW,  IAXIS), &
            ihi => block%limits(HIGH, IAXIS), &
            jlo => block%limits(LOW,  JAXIS), &
            jhi => block%limits(HIGH, JAXIS), &
            klo => block%limits(LOW,  KAXIS), &
            khi => block%limits(HIGH, KAXIS))
    fixedDensity = .FALSE.
    do     k = klo, khi 
      do   j = jlo, jhi 
        do i = ilo, ihi 
          if (solnData(i, j, k, DENS_VAR) < gr_smallrho) then
            solnData(i, j, k, DENS_VAR) = gr_smallrho
            fixedDensity = .TRUE.
          end if
        end do
      end do
    end do
  end associate

  if (fixedDensity) then
    write(*,*) "WARNING: Unphysically small density values set to smallrho"
  end if

  do var = UNK_VARS_BEGIN, UNK_VARS_END
    if (gr_vartypes(var) == VARTYPE_PER_MASS) then
      solnData(:,:,:,var) = solnData(:,:,:,DENS_VAR)*solnData(:,:,:,var)
    end if
  end do
#endif
  
  call Grid_releaseBlkPtr(block, solnData, CENTER)

end subroutine gr_primitiveToConserve

