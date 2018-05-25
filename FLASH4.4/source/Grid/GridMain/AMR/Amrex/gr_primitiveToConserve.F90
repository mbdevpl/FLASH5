!!****if* source/Grid/GridMain/paramesh/gr_primitiveToConserve
!!
!! NAME
!!
!!  gr_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  gr_primitiveToConserve(block_metadata_t(IN) :: block
!!
!!
!! DESCRIPTION
!!
!!  Given a block of data, convert variables which are normally represented in
!!  PER_MASS form (e.g., velocity) to the corresponding conservative form
!!  (i.e., momentum) if gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!
!! ARGUMENTS
!! 
!!   block - the metadata representation of block whose data shall be
!!           transformed
!!
!! NOTES
!!
!!  The variables that are converted are the named cell-centered
!!  solution variables marked to be of type PER_MASS explicitly in a
!!  Config file.  Additionally, abundances and mass scalars are
!!  considered to be of type PER_MASS.
!!
!!  For proper functioning, DENS_VAR must not be marked as PER_MASS!
!!
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_conserveToPrimitive
!!
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are now conserved.  No
!!  such mechanism exists in the code yet.
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

!!REORDER(4):solnData

#include "Flash.h"
#include "constants.h"

subroutine gr_primitiveToConserve(block)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,        ONLY : gr_convertToConsvdForMeshCalls, &
                               gr_vartypes
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN) :: block

  real, pointer :: solnData(:,:,:,:) => null()
  integer       :: i, j, k, var

  if (.NOT. ANY(gr_vartypes == VARTYPE_PER_MASS))    RETURN
  if (.NOT. gr_convertToConsvdForMeshCalls)          RETURN

#ifdef DENS_VAR
  if (gr_vartypes(DENS_VAR) == VARTYPE_PER_MASS) then
    call Driver_abortFlash('[gr_primitiveToConserve] density is PER_MASS')
  end if

  call Grid_getBlkPtr(block, solnData, CENTER)

  ! DEV: TODO Add non-permanent GC here
  associate(ilo => block%limits(LOW,  IAXIS), &
            ihi => block%limits(HIGH, IAXIS), &
            jlo => block%limits(LOW,  JAXIS), &
            jhi => block%limits(HIGH, JAXIS), &
            klo => block%limits(LOW,  KAXIS), &
            khi => block%limits(HIGH, KAXIS))
  do     k = klo, khi 
    do   j = jlo, jhi 
      do i = ilo, ihi 
        if (solnData(DENS_VAR,i,j,k) == 0.0) then
          do var = UNK_VARS_BEGIN, UNK_VARS_END
            if (      (gr_vartypes(var) == VARTYPE_PER_MASS) &
                .AND. (solnData(var,i,j,k) /= 0.0)) then
              ! DEV: TODO Add comment about checking this
              call Driver_abortFlash("[gr_primitiveToConserve] 0/0 error")
            end if
          end do
        end if
      end do
    end do
  end do
  end associate

  do var = UNK_VARS_BEGIN, UNK_VARS_END
    if (gr_vartypes(var) == VARTYPE_PER_MASS) then
      solnData(var,:,:,:) = solnData(DENS_VAR,:,:,:)*solnData(var,:,:,:)
    end if
  end do
 
  call Grid_releaseBlkPtr(block, solnData, CENTER)
#endif

end subroutine gr_primitiveToConserve

