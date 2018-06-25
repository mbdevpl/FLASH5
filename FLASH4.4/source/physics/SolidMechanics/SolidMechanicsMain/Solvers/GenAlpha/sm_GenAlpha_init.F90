#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_GenAlpha_init(ibd,restart)

  use SolidMechanics_data, only : sm_BodyInfo, sm_structure, sm_NumBodies
  use sm_GenAlpha_data, only: sm_GenAlpha_type, sm_GenAlpha_info
  use sm_pk_interface, only: sm_pk_apply, sm_pk_apply_entireBody
  use sm_assemble_interface, only: sm_assemble_mass, sm_assemble_stiff
  use sm_Misc_interface, only: sm_c2fortran_dgssv
  use sm_integinterface, only: sm_ga_initdamp, sm_ga_compute_qddn0, sm_GenAlpha_readCheckpoint
  use Driver_interface, only : Driver_abortFlash, Driver_getSimTime, Driver_getDT
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  
  ! IO Variables
  integer, intent(in) :: ibd
  logical, intent(in) :: restart

  ! Internal Variables
  type(sm_structure), pointer :: body
  type(sm_GenAlpha_type), pointer :: integ
  real :: time, dt
  integer :: neq, ndofs

  ! set Body:
  body  => sm_BodyInfo(ibd)

  ! Check if sm_PredCorr_info has been "allocated"
  if( .not. associated( sm_GenAlpha_info ) ) then
     allocate( sm_GenAlpha_info( sm_NumBodies ) )
  end if

  ! set integ pointer:
  integ => sm_GenAlpha_info(ibd)

  call Driver_getSimTime(time)
  call Driver_getDT(dt)

  !Set their pcconvsflag to TRUE
  integ%pcconvflag = SM_PCNOTCONVERGED

  ! Get the values of epsilon and rhoinf
  call RuntimeParameters_get("pcepsilon",integ%pcepsilon)
  call RuntimeParameters_get("garhoinf",integ%rhoinf)
  call RuntimeParameters_get("gapredepsilon", integ%pred_epsilon)
  call RuntimeParameters_get("gacorrepsilon", integ%corr_epsilon)
  call RuntimeParameters_get("gapredrelax",integ%pred_pos_relax)

  call RuntimeParameters_get("gapredictoronly",integ%predictor_only)
  if( integ%predictor_only .eq. SM_TRUE ) then
     write(*,*) 'GenAlpha: *** WARNING only using predictor ***'
  endif

  call RuntimeParameters_get("gazeroqnonly",integ%zeroQn_only)
  if( integ%zeroQn_only .eq. SM_TRUE ) then
     write(*,*) 'GenAlpha: *** WARNING forcing deformation to be always ZERO ***'
  endif

  ! compute the parameters of the Generalized Alpha Method
  integ%alpha_m = (2*integ%rhoinf - 1.)/(integ%rhoinf + 1.)
  integ%alpha_f = integ%rhoinf/(integ%rhoinf + 1.)
  integ%gamma   = 0.5 - integ%alpha_m + integ%alpha_f
  integ%beta    = 0.25*(1. - integ%alpha_m + integ%alpha_f)**2

  ! Init the needed containers for previous values of Qn
  neq = body%neq
  ndofs = body%ndofs
  if( .not. allocated( body%qi ) ) then
     allocate( body%qi(ndofs) )
     body%qi(1:ndofs) = body%qn(1:ndofs)
  end if
  if( .not. allocated( body%qdi ) ) then
     allocate( body%qdi(ndofs) )
     body%qdi(1:ndofs) = body%qdn(1:ndofs)
  end if
  if( .not. allocated( body%qddi ) ) then
     allocate( body%qddi(ndofs) )
     body%qddi(1:ndofs) = body%qddn(1:ndofs)
  end if
  if( .not. allocated( body%qms ) ) then
     allocate( body%qms(ndofs,1) )
     body%qms = 0.
  end if
  if( .not. allocated( body%Hsn ) ) then
     allocate( body%Hsn(neq) )
     body%Hsn(1:neq) = 0.
  end if
  if( .not. allocated( body%Qsn ) ) then
     allocate( body%Qsn(neq) )
     body%Qsn(1:neq) = 0.
  end if

  ! If this is a restart file, load from checkpoint file
  ! sm_ga_chkpt.#(ibd).#(chkpt).h5
  if( restart ) then

     call sm_GenAlpha_readCheckpoint(ibd)
     
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

  if ( neq .gt. 0 ) then

     ! Mass Matrix for {q,v}
     call sm_assemble_mass(ibd)

     ! Build Damping Matrix (in reference config)
     call sm_ga_initdamp(ibd)
     
     if( .not. restart ) then

        ! Internal solid force vector Fint and Stiffness Matrix K
        call sm_assemble_stiff(ibd, SM_IOPT_QN)

        ! Compute qddn(0)
        call sm_ga_compute_qddn0(ibd)

        ! Copy values of eqn to (n-1)
        body%qddi(1:ndofs) = body%qddn(1:ndofs)
        body%Hsn(1:neq)  = body%Hs(1:neq)
        body%Qsn(1:neq)  = body%Qs(1:neq)

     end if

  endif

  ! Initialize the predictor flag
  integ%pcflag = SM_PCPREDICTOR

  ! Initialize the DT
  integ%dt = dt  

end subroutine sm_GenAlpha_init

