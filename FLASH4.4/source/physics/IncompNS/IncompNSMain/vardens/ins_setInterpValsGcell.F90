!!****if* source/physics/IncompNS/IncompNSMain/vardens/ins_setInterpValsGcell
!!
!!
!! NAME
!!
!!  ins_setInterpValsGcell
!!
!!
!! SYNOPSIS
!!
!!  ins_setInterpValsGcell(logical(IN) :: setval)
!!
!!
!! DESCRIPTION
!!
!! Sets interpolation values for guardcell filling within INS fractional
!! step method.
!!
!!
!! ARGUMENTS
!!
!! setval = .true.  set values 
!!          .false. restore all interpolation-restriction values to quadratic
!!
!!***

subroutine ins_setInterpValsGcell(setval)

#include "Flash.h"

  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work                           
#endif    

  implicit none

  logical, intent(IN) :: setval

  integer, parameter :: intval      = 2 ! Quadratic interpolation and restriction
  integer, parameter :: intfaceval = 32 ! Quadratic restriction + linear restriction
                                        ! in block face.

!  ------------------------------------------------------------------------

#ifdef FLASH_GRID_PARAMESH

  if (setval) then ! Set Gcell interpolation and restriction values:

     ! Work interpolation:
     interp_mask_work(:)= intval;

     ! Cell centered valiables interpolation and Restriction:
     interp_mask_unk(:)     = intval;   
     interp_mask_unk_res(:) = intval;
  
     ! Face Centered variables interpolation and Restriction:
     interp_mask_facex(:) = intval; interp_mask_facex_res(:) = intval;
     interp_mask_facey(:) = intval; interp_mask_facey_res(:) = intval;
     interp_mask_facez(:) = intval; interp_mask_facez_res(:) = intval;

     ! For velocities linear restriction on the face:
     interp_mask_facex_res(VELC_FACE_VAR) = intfaceval
     interp_mask_facey_res(VELC_FACE_VAR) = intfaceval
     interp_mask_facez_res(VELC_FACE_VAR) = intfaceval

  else ! Restore velocity restriction values to quadratic:

     ! For velocities linear restriction on the face:
     interp_mask_facex_res(VELC_FACE_VAR) = intval
     interp_mask_facey_res(VELC_FACE_VAR) = intval
     interp_mask_facez_res(VELC_FACE_VAR) = intval

  endif
  
#endif

  return

end subroutine ins_setInterpValsGcell

