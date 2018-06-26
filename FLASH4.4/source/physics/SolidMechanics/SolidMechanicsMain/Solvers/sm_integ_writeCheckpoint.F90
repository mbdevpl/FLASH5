!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_integ_writeCheckpoint
!!
!!
!! NAME
!!
!! 
!!
!!
!! SYNOPSIS
!!
!!  
!!
!!
!! DESCRIPTION
!! Advances one subiteration, or timestep depending on the time integration desired.
!! For Predictor Corrector-schemes does one predictor or one of the correction 
!! subiterations.
!!
!!***
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_integ_writeCheckPoint()

  use Driver_interface,    only : Driver_abortFlash
  use SolidMechanics_data, only : sm_MeshMe,sm_NumBodies,sm_BodyInfo
  use sm_integinterface,   only : sm_PredCorr_writeCheckpoint, &
                                  sm_Verlet_writeCheckpoint,   &
                                  sm_GenAlpha_writeCheckpoint
  use IO_data,             only: io_checkpointFileNumber

  implicit none

  ! Local Variables
  integer :: ibd
  integer :: sm_checkpt_num

  ! since we are calling the routine BEFORE the number is inc. inside IO_writeCheckpoint
  sm_checkpt_num = io_checkpointFileNumber

  ! Loop over bodies:
  do ibd=1,sm_NumBodies 
     
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
               
        select case(sm_BodyInfo(ibd)%IntegMethod)
        ! Generalized Alpha method
        case(SOLIDINTEG_GENALPHA)
           call sm_GenAlpha_writeCheckpoint(ibd, sm_checkpt_num)
    
        ! Predictor-Corrector method    
        case(SOLIDINTEG_PREDCORR)
           call sm_PredCorr_writeCheckpoint(ibd, sm_checkpt_num)
   
        ! VERLET
        case(SOLIDINTEG_MODVVERLET)
           !call sm_Verlet_writeCheckpoint(ibd,sm_checkpt_num)
            
        case default
           call Driver_abortFlash('SolidMechanics: unknown checkpoint write type.')

        end select

     endif

  enddo ! end looping over ibd


  return

end subroutine sm_integ_writeCheckPoint
