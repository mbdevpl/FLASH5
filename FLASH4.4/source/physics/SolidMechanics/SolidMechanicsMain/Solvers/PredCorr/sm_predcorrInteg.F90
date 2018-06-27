! SUBROUTINE sm_predcorrInteg:
!
!
!======================================================================

subroutine sm_predcorrInteg(ibd,dt)

  use sm_PredCorr_data, only: sm_PredCorr_type, sm_PredCorr_info

  use SolidMechanics_data, only : sm_structure,sm_BodyInfo 

  use sm_integinterface, only : sm_pcem, sm_pcmem, sm_pc_AB2, sm_pc_AM2,   &
                                sm_pc_AB3, sm_pc_AM3, sm_pc_AB4, sm_pc_AM4

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "sm_integrator.h"

  integer, intent(in) :: ibd
  real, intent(in) :: dt

  real, parameter :: CONST_ONE = 1.

  real :: e1
  integer :: neq

  type(sm_structure), pointer :: BodyInfo
  type(sm_PredCorr_type), pointer :: integ

  ! Assign Pointers
  BodyInfo => sm_BodyInfo(ibd)
  integ => sm_PredCorr_info(ibd)
  
  ! number of unrestrained dofs
  neq = BodyInfo % neq

  if (neq .eq. 0) return ! Only unrestrained dofs are time-advanced here. 

  select case(integ%pciter) ! EULER = 1, AB2 = 2, AB3 = 3, HMG = 4

  case(SM_PCEULER)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then
        ! Y,YDOT(-1)                         Y, YDOT (0)
        BodyInfo %  qms(1:neq,-1) = BodyInfo  %qn(1:neq) ! Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-1) = BodyInfo %qdn(1:neq)
        BodyInfo %qddms(1:neq,-1) = BodyInfo%qddn(1:neq)

        ! shift the dt's
        integ%vardt(-1) = integ%vardt(0)
        integ%vardt( 0) = dt


! ...   Predict the Solution ( Euler Method )
!                        Y0            F0                  pY1
        ! Positions q:
        call sm_pcem(neq,BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1),dt,BodyInfo%qn(1:neq))
        ! Velocities qd:
        call sm_pcem(neq,BodyInfo%qdms(1:neq,-1),BodyInfo%qddms(1:neq,-1),dt,BodyInfo%qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

! ...   Correct the Solution ( Modified Euler Method )
!                         Y0               F0              kF1
        ! Positions q:
        call sm_pcmem (neq,BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1),BodyInfo%qdn(1:neq), &
                     dt, BodyInfo%qi(1:neq))
!                                   k+1Y1       
        ! Velocities qd:
        call sm_pcmem (neq,BodyInfo%qdms(1:neq,-1),BodyInfo%qddms(1:neq,-1),BodyInfo%qddn(1:neq), &
                     dt, BodyInfo%qdi(1:neq))
!                                   k+1Y1                  
!       Check for convergence
!
!       e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(BodyInfo %qi(1:neq)-BodyInfo %qn(1:neq))), &
                 maxval(abs(BodyInfo%qdi(1:neq)-BodyInfo%qdn(1:neq)))) 
        

!       kY1            k+1Y1  
        BodyInfo% qn(1:neq) = BodyInfo %qi(1:neq)
        BodyInfo%qdn(1:neq) = BodyInfo%qdi(1:neq)

!       Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCEULER) integ%pciter = SM_PCAB2
        endif

     endif
        
  case(SM_PCAB2)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        BodyInfo %  qms(1:neq,-2) = BodyInfo  %qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-2) = BodyInfo %qdms(1:neq,-1)
        BodyInfo %qddms(1:neq,-2) = BodyInfo%qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        BodyInfo %  qms(1:neq,-1) = BodyInfo  %qn(1:neq) ! Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-1) = BodyInfo %qdn(1:neq)
        BodyInfo %qddms(1:neq,-1) = BodyInfo%qddn(1:neq)

        ! shift the dt's
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB2 )
        ! Positions q:
        call sm_pc_AB2(neq,BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1),   &
                       BodyInfo%qdms(1:neq,-2),integ%vardt(-1:0),1,BodyInfo%qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB2(neq,BodyInfo%qdms(1:neq,-1),BodyInfo%qddms(1:neq,-1), &
                       BodyInfo%qddms(1:neq,-2),integ%vardt(-1:0),1,BodyInfo%qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 2 )
        ! Positions q:
        call sm_pc_AM2(neq,BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1), &
                           BodyInfo%qdms(1:neq,-2),BodyInfo%qdn(1:neq),    &
                           integ%vardt(-1:0),1, BodyInfo%qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM2 (neq, BodyInfo%qdms(1:neq,-1),  BodyInfo%qddms(1:neq,-1), &
                             BodyInfo%qddms(1:neq,-2), BodyInfo%qddn(1:neq), &
                             integ%vardt(-1:0),1, BodyInfo%qdi(1:neq))

        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(BodyInfo %qi(1:neq)-BodyInfo %qn(1:neq))), &
                 maxval(abs(BodyInfo%qdi(1:neq)-BodyInfo%qdn(1:neq)))) 
        

        ! kY1            k+1Y1  
        BodyInfo% qn(1:neq) = BodyInfo %qi(1:neq)
        BodyInfo%qdn(1:neq) = BodyInfo%qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCAB2) integ%pciter = SM_PCAB3
        endif

     endif

  case(SM_PCAB3)

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-3)                         Y, YDOT (-2)
        BodyInfo %  qms(1:neq,-3) = BodyInfo  %qms(1:neq,-2) ! n-3 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-3) = BodyInfo %qdms(1:neq,-2)
        BodyInfo %qddms(1:neq,-3) = BodyInfo%qddms(1:neq,-2)

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        BodyInfo %  qms(1:neq,-2) = BodyInfo  %qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-2) = BodyInfo %qdms(1:neq,-1)
        BodyInfo %qddms(1:neq,-2) = BodyInfo%qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        BodyInfo %  qms(1:neq,-1) = BodyInfo  %qn(1:neq) ! Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-1) = BodyInfo %qdn(1:neq)
        BodyInfo %qddms(1:neq,-1) = BodyInfo%qddn(1:neq)

        ! shift the dt's
        integ%vardt(-3) = integ%vardt(-2)
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB3 )
        ! Positions q:
        call sm_pc_AB3(neq,BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1),   &
                           BodyInfo%qdms(1:neq,-2), BodyInfo%qdms(1:neq,-3), &
                           integ%vardt(-2:0),2,BodyInfo%qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB3(neq, BodyInfo%qdms(1:neq,-1),BodyInfo%qddms(1:neq,-1),  &
                            BodyInfo%qddms(1:neq,-2),BodyInfo%qddms(1:neq,-3), &
                            integ%vardt(-2:0),2,BodyInfo%qdn(1:neq))        

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 3 )
        ! Positions q:
        call sm_pc_AM3(neq, BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1), &
                            BodyInfo%qdms(1:neq,-2), BodyInfo%qdms(1:neq,-3), &
                            BodyInfo%qdn(1:neq),    &
                            integ%vardt(-2:0),2, BodyInfo%qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM3(neq, BodyInfo%qdms(1:neq,-1),  BodyInfo%qddms(1:neq,-1), &
                            BodyInfo%qddms(1:neq,-2), BodyInfo%qddms(1:neq,-3), &
                            BodyInfo%qddn(1:neq), &
                            integ%vardt(-2:0),2, BodyInfo%qdi(1:neq))

        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(BodyInfo %qi(1:neq)-BodyInfo %qn(1:neq))), &
                 maxval(abs(BodyInfo%qdi(1:neq)-BodyInfo%qdn(1:neq)))) 
        
        ! kY1            k+1Y1  
        BodyInfo% qn(1:neq) = BodyInfo %qi(1:neq)
        BodyInfo%qdn(1:neq) = BodyInfo%qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           if (integ%pcmethod .gt. SM_PCAB3) integ%pciter = SM_PCAB4
        endif

     endif

  case(SM_PCAB4)        

     if (integ%pcflag .eq. SM_PCPREDICTOR) then

        ! Y,YDOT(-4)                         Y, YDOT (-3)
        BodyInfo %  qms(1:neq,-4) = BodyInfo  %qms(1:neq,-3) ! n-4 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-4) = BodyInfo %qdms(1:neq,-3)
        BodyInfo %qddms(1:neq,-4) = BodyInfo%qddms(1:neq,-3)

        ! Y,YDOT(-3)                         Y, YDOT (-2)
        BodyInfo %  qms(1:neq,-3) = BodyInfo  %qms(1:neq,-2) ! n-3 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-3) = BodyInfo %qdms(1:neq,-2)
        BodyInfo %qddms(1:neq,-3) = BodyInfo%qddms(1:neq,-2)

        ! Y,YDOT(-2)                         Y, YDOT (-1)
        BodyInfo %  qms(1:neq,-2) = BodyInfo  %qms(1:neq,-1) ! n-2 Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-2) = BodyInfo %qdms(1:neq,-1)
        BodyInfo %qddms(1:neq,-2) = BodyInfo%qddms(1:neq,-1)

        ! Y,YDOT(-1)                         Y, YDOT (0)
        BodyInfo %  qms(1:neq,-1) = BodyInfo  %qn(1:neq) ! Previous steps pos, vel, accel
        BodyInfo % qdms(1:neq,-1) = BodyInfo %qdn(1:neq)
        BodyInfo %qddms(1:neq,-1) = BodyInfo%qddn(1:neq)

        integ%vardt(-4) = integ%vardt(-3)
        integ%vardt(-3) = integ%vardt(-2)
        integ%vardt(-2) = integ%vardt(-1)
        integ%vardt(-1) = integ%vardt( 0)
        integ%vardt( 0) = dt

        ! ...   Predict the Solution ( AB4 )
        ! Positions q:
        call sm_pc_AB4(neq, BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1),   &
                            BodyInfo%qdms(1:neq,-2), BodyInfo%qdms(1:neq,-3), &
                            BodyInfo%qdms(1:neq,-4), &
                            integ%vardt(-3:0),3,BodyInfo%qn(1:neq))
        ! Velocities qd:
        call sm_pc_AB4(neq, BodyInfo%qdms(1:neq,-1),BodyInfo%qddms(1:neq,-1),  &
                            BodyInfo%qddms(1:neq,-2),BodyInfo%qddms(1:neq,-3), &
                            BodyInfo%qddms(1:neq,-4), &
                            integ%vardt(-3:0),3,BodyInfo%qdn(1:neq))  

        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = CONST_ONE !!! We are not ok with an error of 1. !!!
        
     else

        ! ...   Correct the Solution ( Adams-Moulton 4 )
        ! Positions q:
        call sm_pc_AM4(neq, BodyInfo%qms(1:neq,-1),BodyInfo%qdms(1:neq,-1), &
                            BodyInfo%qdms(1:neq,-2), BodyInfo%qdms(1:neq,-3), &
                            BodyInfo%qdms(1:neq,-4), &
                            BodyInfo%qdn(1:neq),     &
                            integ%vardt(-3:0),3, BodyInfo%qi(1:neq))
      
        ! Velocities qd:
        call sm_pc_AM4(neq, BodyInfo%qdms(1:neq,-1),  BodyInfo%qddms(1:neq,-1), &
                            BodyInfo%qddms(1:neq,-2), BodyInfo%qddms(1:neq,-3), &
                            BodyInfo%qddms(1:neq,-4), &
                            BodyInfo%qddn(1:neq), &
                            integ%vardt(-3:0),3, BodyInfo%qdi(1:neq))
        
        ! Check for convergence
        ! e1 = L1 norm of k+1Y1-kY1 
        e1 = max(maxval(abs(BodyInfo %qi(1:neq)-BodyInfo %qn(1:neq))), &
                 maxval(abs(BodyInfo%qdi(1:neq)-BodyInfo%qdn(1:neq)))) 
        
        ! kY1            k+1Y1  
        BodyInfo% qn(1:neq) = BodyInfo %qi(1:neq)
        BodyInfo%qdn(1:neq) = BodyInfo%qdi(1:neq)

        ! Load convergence error for the body:        
        integ%pcerr = e1

        if (e1 .le. integ%pcepsilon) then
           integ%pcflag      = SM_PCPREDICTOR
           integ%pcconvflag  = SM_PCCONVERGED
           !!if (integ%pcmethod .gt. SM_PCAB4) integ%pciter = SM_PCAB5
        endif

     endif
        

  case default
     
     call Driver_abortFlash("Error in sm_predcorrInteg: sm_pciter not equal 1-4.")
        
  end select

  nullify(BodyInfo)
  nullify(integ)

  return

end subroutine sm_predcorrInteg
