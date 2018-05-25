!!****if* source/Grid/GridMain/paramesh/gr_conserveToPrimitive
!!
!! NAME
!!
!!  gr_conserveToPrimitive
!!
!!
!! SYNOPSIS
!!
!!  gr_conserveToPrimitive(block_metadata_t(in) :: block,
!!                         logical(in)          :: allCells)
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a block of data, convert variables which are normally
!!  represented in PER_MASS form (e.g., velocity) from the 
!!  corresponding conservative form (i.e., momentum) back to the normal
!!  PER_MASS form  if gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!  Additionally,
!!   - energies (ENER_VAR and EINT_VAR) are forced to be .ge. gr_smalle,
!!   - the density (DENS_VAR) is forced to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!! ARGUMENTS
!! 
!!   block - the metadata representation of block whose data shall be
!!           transformed
!!
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
!!
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
!!
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_primitiveToConserve 
!!
!!
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are no longer conserved.  No
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
        if (needToConvert) then
          if (solnData(DENS_VAR,i,j,k) == 0.0) then
            dens_old_inv = 0.0
            do var = UNK_VARS_BEGIN, UNK_VARS_END
              if (      (gr_vartypes(var) == VARTYPE_PER_MASS) &
                  .AND. (solnData(var,i,j,k) /= 0.0)) then
                write(*,*) "0/0 GC error", i, j, k, solnData(var,i,j,k)
                call Driver_abortFlash("[gr_conserveToPrimitive] 0/0 GC Error")
              end if
            end do
          else
            dens_old_inv = 1.0 / solnData(DENS_VAR,i,j,k)
          end if

          do var = UNK_VARS_BEGIN, UNK_VARS_END
            if (gr_vartypes(var) == VARTYPE_PER_MASS) then
              solnData(var,i,j,k) = dens_old_inv * solnData(var,i,j,k)
            end if
          end do
        end if

        ! small limits -- in case the interpolants are not monotonic
        solnData(DENS_VAR,i,j,k) = max(solnData(DENS_VAR,i,j,k), gr_smallrho)
      end do
    end do
  end do
#endif

#ifdef ENER_VAR               
  solnData(ENER_VAR,:,:,:) = max(solnData(ENER_VAR,:,:,:), gr_smalle)
#endif
#ifdef EINT_VAR
  solnData(EINT_VAR,:,:,:) = max(solnData(EINT_VAR,:,:,:), gr_smalle)
#endif

  call Grid_releaseBlkPtr(block, solnData, CENTER)

end subroutine gr_conserveToPrimitive
        
