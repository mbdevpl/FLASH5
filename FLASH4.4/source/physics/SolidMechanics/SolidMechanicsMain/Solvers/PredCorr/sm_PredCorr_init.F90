#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_PredCorr_init(ibd,restart)

      use SolidMechanics_data, only : sm_BodyInfo, sm_structure, sm_NumBodies
      use sm_PredCorr_data, only: sm_PredCorr_type, sm_PredCorr_info
      use RuntimeParameters_interface, ONLY : RuntimeParameters_get
      use sm_pk_interface, only: sm_pk_apply, sm_pk_apply_entireBody, sm_pk_updatekinematics_rigid
      use Driver_interface, only : Driver_abortFlash, Driver_getSimTime, Driver_getDT
      use sm_integinterface, only : sm_pc_initmass, sm_PredCorr_readCheckpoint
      use sm_assemble_interface, only: sm_assemble_IntForce

      implicit none

      ! IO Variables
      integer, intent(in) :: ibd
      logical, intent(in) :: restart

      ! Internal Variables
      type(sm_structure), pointer :: body
      type(sm_PredCorr_type), pointer :: integ
      real :: time, dt
      integer :: neq

      ! set Body:
      body  => sm_BodyInfo(ibd)

      ! Check if sm_PredCorr_info has been "allocated"
      if( .not. associated( sm_PredCorr_info ) ) then
         allocate( sm_PredCorr_info( sm_NumBodies ) )
      end if

      ! set integ pointer:
      integ => sm_PredCorr_info(ibd)

      call Driver_getSimTime(time)
      call Driver_getDT(dt)

      !Set their pcconvsflag to TRUE
      integ%pcconvflag = SM_PCNOTCONVERGED
      call RuntimeParameters_get("pcmethod",integ%pcmethod)
      call RuntimeParameters_get("pcepsilon",integ%pcepsilon)

      ! Init the needed containers for previous values of Qn
      neq = body%neq
      allocate( body%  qms(1:neq, -integ%pcmethod:-1), &
                body% qdms(1:neq, -integ%pcmethod:-1), &
                body%qddms(1:neq, -integ%pcmethod:-1), &
                body%qdi(1:neq) )

      ! Init the dt tracker for variable timesteps
      allocate( integ%vardt(-integ%pcmethod:0) )
      integ%vardt(:) = 0.
      integ%vardt(0) = dt

      ! If this is a restart file, load from checkpoint file
      ! sm_ga_chkpt.#(ibd).#(chkpt).h5
      if( restart ) then

         call sm_PredCorr_readCheckpoint(ibd)

      else
         !Specified Coordinates
         ! This updates the values in restrained part of qn,qdn,qddn
         if( body%ng .gt. 0 ) then ! can be replaced by another flag if needed
         
            if( body%IC_flag == SM_IC_PK ) then
               call sm_pk_apply_entireBody(ibd, time)
            else
               call sm_pk_apply(ibd,time)
            end if

         end if

      end if 

      ! Initialize the sm_pciter
      if(.not. restart) integ%pciter = SM_PCEULER

      ! Initialize the predictor flag
      integ%pcflag = SM_PCPREDICTOR

      ! Initialize the DT
      integ%dt = dt

      ! Set initial configuration for rigid bodies:
      if( body%bodyType .eq. BODYTYPE_RIGID)  then
        call sm_pk_updatekinematics_rigid(ibd,restart,SM_INIT)
        return
      endif  

      if ( neq .gt. 0 ) then
        ! Mass Matrix for qi,qv
        call sm_pc_initmass(ibd)

        ! Build Damping Matrix (in reference config)
        call sm_pc_initdamp(ibd)
          
        ! Internal solid force vector Fint 
        call sm_assemble_IntForce(ibd)
      endif



end subroutine sm_PredCorr_init

