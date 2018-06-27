!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_integ_advance
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

subroutine sm_integ_advance(restart_local)

  use Driver_interface, only : Driver_abortFlash
  use SolidMechanics_data, only : sm_MeshMe,sm_NumBodies,sm_BodyInfo, sm_structure
  use sm_integinterface, only : sm_PredCorr_advance, sm_Verlet_advance, sm_GenAlpha_advance
  use sm_contact_interface, only : sm_contact

  implicit none

#include "Flash_mpi.h"

  ! IO Variables
  logical, INTENT(IN) :: restart_local

  ! Local Variables
  integer :: ibd

  ! Meta-body gather buffer
  integer :: ierr, ndofs, MetaBody0

  ! Loop over bodies:
  do ibd=1,sm_NumBodies

     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

        select case(sm_BodyInfo(ibd)%IntegMethod)

           ! Generalized Alpha method
        case(SOLIDINTEG_GENALPHA)
           call sm_GenAlpha_advance(ibd,restart_local)

           ! Predictor-Corrector method
        case(SOLIDINTEG_PREDCORR)
           call sm_PredCorr_advance(ibd,restart_local)

           ! VERLET
        case(SOLIDINTEG_MODVVERLET)
           call sm_Verlet_advance(ibd,restart_local)

        case default
           call Driver_abortFlash('SolidMechanics: Unknown advance integrtion scheme.')

        end select

        !
        ! sync meta-bodies, call an MPI_scatter on metabody 0
        !
        if ( sm_BodyInfo(ibd)%MetaBody .gt. CONSTANT_ZERO ) then

           ndofs = sm_BodyInfo(ibd)%ndofs
           MetaBody0 = 0

           ! qi:
           call MPI_Bcast( sm_BodyInfo(ibd)%qi, ndofs, FLASH_REAL, MetaBody0, &
                sm_BodyInfo(ibd)%mbcomm, ierr)

        end if

     end if


  enddo ! end looping over ibd


  if (sm_Numbodies>1) then
     do ibd=1,sm_NumBodies
        call sm_contact(ibd, restart_local)
     end do

  end if

  return

end subroutine sm_integ_advance
