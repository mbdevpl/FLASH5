!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_ioInit
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
!!  Initialize IO variables and files for SolidMechanics.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_ioInit(restart)

  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo
  use sm_iointerface, only : sm_ioInit_rigid,sm_ioInit_3DFlexible,sm_ioInit_RBC
  use Driver_data, only : dr_simTime

  implicit none

  ! IO Variables
  logical, INTENT(IN) :: restart

  ! Local Vars:
  integer :: ibd

  do ibd=1,sm_NumBodies
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

        ! Initialize IO based on BodyType
        select case( sm_BodyInfo(ibd)%BodyType )

        case( BODYTYPE_RIGID )         
           call sm_ioInit_rigid(restart,ibd,dr_simTime)

#if NDIM == MDIM

        case( BODYTYPE_3DFLEXIBLE )
           call sm_ioInit_3DFlexible(restart,ibd,dr_simTime)

        case( BODYTYPE_RBC )
           call sm_ioInit_RBC(restart,ibd,dr_simTime)

#else /* 2D */

        case( BODYTYPE_2DFLEXIBLE )
           call Driver_abortFlash('2D FLEXIBLE BodyType IO not yet implemented')

#endif

        case default
           call Driver_abortFlash("sm_ioInit : BodyType IO not yet implemented")

        end select

     endif
  enddo

  return

end subroutine sm_ioInit
