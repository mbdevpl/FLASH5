!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_ioWriteStates
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
!!  Write variables of interest from SolidMechanics.
!!  This routine will append the current vars of interest
!!  To the files initialized in sm_ioInit.F90.
!!
!!***

subroutine sm_ioWriteStates()

  use SolidMechanics_data, only : sm_MeshMe, sm_NumBodies, sm_BodyInfo
  use sm_iointerface, only : sm_ioWriteStates_rigid,       &
                             sm_ioWriteStates_3DFlexible,  &
                             sm_ioWriteStates_rbc
  use Driver_data, only : dr_simTime, dr_nstep

  implicit none
#include "SolidMechanics.h"

  ! Local Vars:
  integer :: ibd

  do ibd=1,sm_NumBodies
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

        ! Initialize IO based on BodyType
        select case( sm_BodyInfo(ibd)%BodyType )

        case( BODYTYPE_RIGID )         
           call sm_ioWriteStates_rigid(ibd,dr_nstep,dr_simTime)

        case( BODYTYPE_2DFLEXIBLE )
           call Driver_abortFlash('2D FLEXIBLE BodyType IO not yet implemented')

        case( BODYTYPE_3DFLEXIBLE )
           call sm_ioWriteStates_3DFlexible(ibd,dr_nstep,dr_simTime)

        case( BODYTYPE_RBC )
           call sm_ioWriteStates_rbc(ibd,dr_nstep,dr_simTime)

        case default
           call Driver_abortFlash("sm_ioWriteStates : BodyType IO not yet implemented")

        end select

     endif
  enddo

  return


end subroutine sm_ioWriteStates

