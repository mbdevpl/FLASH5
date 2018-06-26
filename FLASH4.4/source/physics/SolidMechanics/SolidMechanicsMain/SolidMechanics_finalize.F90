!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SolidMechanics_finalize
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine SolidMechanics_finalize()

  USE HDF5
  use sm_pk_interface, only: sm_pk_finalize
  use sm_Misc_interface, only: sm_deallocateBody
  use SolidMechanics_data, only : sm_MeshMe, sm_NumBodies, sm_bodyInfo
  use Driver_interface, only : Driver_getMype

  implicit none

#include "constants.h"
#include "SolidMechanics.h"

  integer             :: h5ferr, ibd

  ! Which Processor am I:
  call Driver_getMype(GLOBAL_COMM, sm_MeshMe)
  

  if (sm_MeshMe .eq. MASTER_PE) write(*,*)'deallocating sm_bodyInfo...'
  ! close all body related stuff
  do ibd=1,sm_NumBodies
     if (sm_MeshMe== sm_bodyInfo(ibd)%BodyMaster ) then 
        ! kill each body
        call sm_deallocateBody(ibd)
        
        ! close the Prescribed kinematics stuff
        call sm_pk_finalize()
     end if
  end do
  if( associated(sm_bodyInfo) ) deallocate( sm_bodyInfo )

  if (sm_MeshMe .eq. MASTER_PE) write(*,*)'Done.'
  
  ! close the HDF5 fortran library
  call h5close_f(h5ferr)

  return

end subroutine SolidMechanics_finalize
