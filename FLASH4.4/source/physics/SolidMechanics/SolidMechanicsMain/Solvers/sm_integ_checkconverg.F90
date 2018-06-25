!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/sm_integ_checkconverg
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
!! Checks for convergence on schemes which iterate to convergence on a time step.
!!
!!***
#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"
#include "sm_integrator.h"


subroutine sm_integ_checkconverg(convflag_all)

  use SolidMechanics_data, only : sm_MeshMe,sm_NumBodies,sm_BodyInfo
  use sm_integdata, only : sm_convflag_all, sm_errmax_all, sm_integ_subiter, sm_err_diverge, sm_integ_subiter_old
  use sm_integinterface, only : sm_PredCorr_checkconverg, sm_GenAlpha_checkconverg

  implicit none
  include "Flash_mpi.h"  
  integer, INTENT(OUT) :: convflag_all

  ! Local Variables   
  integer :: convflag_loc, convflag_pc_loc
  integer :: ibd, ierr
  real :: errmax_loc, errmax_pc_loc

  convflag_loc = SM_CONVERGED
  errmax_loc   = 0.
  convflag_all = SM_NOTCONVERGED
  
  ! First test if any proc has a body wnich in not converged 
  ! (convflag_loc == SM_NOTCONVERGED)
  do ibd=1,sm_NumBodies
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

        select case(sm_BodyInfo(ibd)%IntegMethod)
        ! Generalized Alpha method
        case(SOLIDINTEG_GENALPHA)

           if (sm_BodyInfo(ibd)%neq .gt. 0) then ! If there are unknowns check convergence

              call sm_GenAlpha_checkconverg(SM_TESTCONVERGE,ibd,convflag_pc_loc,errmax_pc_loc)
              convflag_loc  = convflag_pc_loc * convflag_loc
              errmax_loc    = max(errmax_loc , errmax_pc_loc)

           endif
   
        ! Adams Predictor-Corrector method    
        case(SOLIDINTEG_PREDCORR)

           if (sm_BodyInfo(ibd)%neq .gt. 0) then ! If there are unknowns check convergence

              call sm_PredCorr_checkconverg(SM_TESTCONVERGE,ibd,convflag_pc_loc,errmax_pc_loc)
              convflag_loc  = convflag_pc_loc * convflag_loc
              errmax_loc    = max(errmax_loc , errmax_pc_loc)

           endif
  
        ! VERLET
        case(SOLIDINTEG_MODVVERLET)
           !call Driver_abortFlash('sm_integ_checkconverg: Verlet not implemented.')

        case default
           call Driver_abortFlash('sm_integ_checkconverg: Unknown advance integrtion scheme.')

        end select

     endif

  enddo ! end looping over ibd


  ! Do a logical Reduce of convflag_loc, max error 
  call MPI_ALLREDUCE(convflag_loc,sm_convflag_all,CONSTANT_ONE,FLASH_INTEGER, &
                     MPI_PROD,FLASH_COMM,ierr)

  call MPI_ALLREDUCE(errmax_loc,sm_errmax_all,CONSTANT_ONE,FLASH_REAL, &
                     MPI_MAX,FLASH_COMM,ierr)  


  ! If one body is not converged, all bodies perform one more sub-iteration:
  if (sm_convflag_all .eq. SM_NOTCONVERGED) then

     do ibd=1,sm_NumBodies
        if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

           select case(sm_BodyInfo(ibd)%IntegMethod)

           ! Generalized Alpha method
           case(SOLIDINTEG_GENALPHA)
              ! Here convflag_pc_loc,errmax_pc_loc are not used.
              call sm_GenAlpha_checkconverg(SM_SETNOTCONVERGED,ibd,convflag_pc_loc,errmax_pc_loc)

           ! Adams Predictor-Corrector method    
           case(SOLIDINTEG_PREDCORR)

              ! Here convflag_pc_loc,errmax_pc_loc are not used.
              call sm_PredCorr_checkconverg(SM_SETNOTCONVERGED,ibd,convflag_pc_loc,errmax_pc_loc)
                                                                                            
           ! VERLET
           case(SOLIDINTEG_MODVVERLET)
              call Driver_abortFlash('sm_integ_checkconverg: Verlet not implemented.')

           case default
              call Driver_abortFlash('sm_integ_checkconverg: Unknown advance integrtion scheme.')

           end select

        endif

     enddo ! end looping over ibd

     ! Write to screen
     sm_integ_subiter = sm_integ_subiter + 1
     if(sm_MeshMe .eq. MASTER_PE) then
        write(*,'(A)') '    -------------------------------------------'
        write(*,'(A,I3,A,g18.10)') '    SubIt=',sm_integ_subiter,' ErrMax=',sm_errmax_all
     end if
       
     if(sm_errmax_all .ge. sm_err_diverge) then
        call Driver_abortFlash('sm_integ_checkconverg : Sub-iteration diverged.')
     endif


  else ! All bodies converged
 
     sm_integ_subiter_old = sm_integ_subiter
     sm_integ_subiter = 0
     if(sm_MeshMe .eq. MASTER_PE) write(*,*) 'Step Converged. ErrMax=',sm_errmax_all

  endif


  convflag_all = sm_convflag_all

  return

end subroutine sm_integ_checkconverg
