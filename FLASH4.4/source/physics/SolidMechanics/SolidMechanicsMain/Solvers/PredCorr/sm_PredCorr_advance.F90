! SolidMechanicsMain/Solvers/PredCorr

#include "SolidMechanics.h"

subroutine sm_PredCorr_advance(ibd,restart)
  
  use SolidMechanics_data, only: sm_BodyInfo
  use sm_PredCorr_data, only: sm_PredCorr_info
  use sm_integinterface, only: sm_pc_compute_qddn, sm_predcorrInteg
  use sm_pk_interface, only: sm_pk_apply, sm_pk_updatekinematics_rigid, &
                             sm_pk_angvelconstraint_rigid
  use sm_assemble_interface, only: sm_assemble_mass, sm_assemble_IntForce, sm_assemble_ExtForce
  use Driver_interface, only: Driver_getSimTime

  implicit none
#include "sm_integrator.h"
  ! IO Variables
  integer, intent(in) :: ibd
  logical, intent(in) :: restart

  ! Internal Variables
  real :: time

  ! Get the simulation time
  call Driver_getSimTime(time)

  ! Get the current DT
  call Driver_getDT(sm_PredCorr_info(ibd)%dt)

  ! Specified Coordinates
  ! This updates the values in restrained part of qn,qdn,qddn
  if( sm_BodyInfo(ibd)%ng .gt. 0 ) then ! There are restraints
     call sm_pk_apply(ibd,time)
  end if

  ! Fluid Forces + External Body Forces Fext
  ! This populates body%Hs
  call sm_assemble_ExtForce(ibd,time)

  if( sm_BodyInfo(ibd)%neq .eq. 0 ) then
     ! Case of rigid bodies: Populate orientation, and angular velocity-acceleration data:
     if( sm_BodyInfo(ibd)%bodyType .eq. BODYTYPE_RIGID ) &
         call sm_pk_updatekinematics_rigid(ibd,.false.,SM_ADVANCE)
     return ! There are no unknowns
  endif

  ! Mass Matrix for qi,qv
  if(  sm_BodyInfo(ibd)%flag_constMass == SM_FALSE ) then
     call sm_assemble_mass(ibd)
  end if

  ! Internal solid force vector Fint, 
  ! this populates body%Qs(1:neq)
  call sm_assemble_IntForce(ibd, SM_IOPT_QN)

  ! Solve for BodyInfo%qdn , qddn
  ! Here BodyInfo%qdn does not need to be updated, its done on sm_predcorrInteg
  ! Solve for BodyInfo%qddn
  call sm_pc_compute_qddn(ibd)

  ! Do the actual integration step for new  BodyInfo%qn,qdn
  call sm_predcorrInteg(ibd,sm_PredCorr_info(ibd)%dt)

  ! Case of rigid bodies: 
  ! 1. If body free to rotate in 3D: Force value of orientation vars rate of change 
  !    with the resulting angular velocity.
  ! 2. Populate orientation, and angular velocity-acceleration data:
  if( sm_BodyInfo(ibd)%bodyType .eq. BODYTYPE_RIGID ) then
      if (sm_PredCorr_info(ibd)%pcconvflag .eq. SM_PCCONVERGED) call sm_pk_angvelconstraint_rigid(ibd)
      call sm_pk_updatekinematics_rigid(ibd,.false.,SM_ADVANCE)
   endif

end subroutine sm_PredCorr_advance

