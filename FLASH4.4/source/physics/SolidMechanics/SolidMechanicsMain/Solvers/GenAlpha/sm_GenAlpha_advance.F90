!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Solvers/GenAlpha/sm_GenAlpha_advance
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
!!
!!
!!
!!***

#include "SolidMechanics.h"
#include "sm_integrator.h"

subroutine sm_GenAlpha_advance(ibd,restart)
  use SolidMechanics_data,   only: sm_BodyInfo, sm_structure
  use sm_GenAlpha_data,      only: sm_GenAlpha_info, sm_GenAlpha_type
  use sm_pk_interface,       only: sm_pk_apply
  use sm_assemble_interface, only: sm_assemble_mass, sm_assemble_stiff, sm_assemble_ExtForce
  use sm_Misc_interface,     only: sm_c2fortran_dgssv,  DGEMV_sparse
  use Driver_interface,      only: Driver_getSimTime, Driver_getDT, Driver_getNstep
  use sm_integdata,          only: sm_integ_subiter
  use sm_integinterface,     only: sm_ga_modstiff_fluid 

  implicit none

  ! IO Variables
  integer, intent(in) :: ibd
  logical, intent(in) :: restart

  ! Internal Variables
  type(sm_structure), pointer :: body
  type(sm_GenAlpha_type), pointer :: integ
  real :: time, dt, a1, a2, a3, b1, b2, b3
  integer :: neq, ndofs, ng, ng1, nnz
  real :: e1, error_q
  real, allocatable, dimension(:) :: vec
  integer :: iopt, ldb, info, nrhs, nstep

  ! Local subiteration variables
  real    :: error_qd, iter_epsilon
  integer :: Nmax_iter, iter_flag, iter

  ! Setup pointers
  body  => sm_BodyInfo(ibd)
  integ => sm_GenAlpha_info(ibd)
  neq   = body%neq
  ndofs = body%ndofs
  ng    = body%ng
  ng1   = neq+1
  allocate(vec(ndofs))

  ! Get the simulation time
  call Driver_getSimTime(time)

  ! Get the current DT
  call Driver_getDT(dt)
  integ%dt = dt

  ! Save initial states
  body%qms(1:ndofs,1) = body%qn(1:ndofs)

  ! If Predictor
  if (integ%pcflag .eq. SM_PCPREDICTOR) then
     ! Save old values
     body%qi(1:ndofs)   = body%qn(1:ndofs)
     body%qdi(1:ndofs)  = body%qdn(1:ndofs)
     body%qddi(1:ndofs) = body%qddn(1:ndofs)
     ! Old Forces
     body%Qsn(1:neq) = body%Qs(1:neq)
     body%Hsn(1:neq) = body%Hs(1:neq)
     ! Epsilon
     iter_epsilon = integ%pred_epsilon
  else
     ! Epsilon
     iter_epsilon = integ%corr_epsilon
  end if


  ! Specified Coordinates
  ! This updates the values in restrained part of qn,qdn,qddn
  if( sm_BodyInfo(ibd)%ng .gt. 0 ) then ! There are restraints
     call sm_pk_apply(ibd,time)
  end if

  ! Fluid Forces + External Body Forces Fext
  ! This populates body%Hs
  call sm_assemble_ExtForce(ibd,time)


  ! Override with zeros
  if( integ%zeroQn_only == SM_TRUE ) then
     body%qn(1:neq)   = 0.
     body%qdn(1:neq)  = 0.
     body%qddn(1:neq) = 0.
  end if
  
  iter_flag = 0
  iter      = 1
  Nmax_iter = 200
  do while( iter_flag == 0 .and. integ%zeroQn_only == SM_FALSE ) 

     ! Mass Matrix for qi,qv
     if(  sm_BodyInfo(ibd)%flag_constMass == SM_FALSE ) then
        call sm_assemble_mass(ibd)
     end if

     ! Internal solid force vector Fint, 
     ! this populates body%Qs(1:neq) and body%K
     call sm_assemble_stiff(ibd, SM_IOPT_QN)

     ! Modified K for implicit pressure calcs:
     !if( integ%pcflag .eq. SM_PCCORRECTOR ) call sm_ga_modstiff_fluid(ibd)


     !-----------------------------------------------------
     ! Build the Residual of the GenAlpha Equation
     !
     a1 = (1.-integ%alpha_m)/(integ%beta*dt**2)
     a2 = -(1.-integ%alpha_m)/(integ%beta*dt)
     a3 = -(1.-integ%alpha_m-2*integ%beta)/(2.*integ%beta)
     b1 = ( (1.-integ%alpha_f)*integ%gamma )/( integ%beta * dt)
     b2 = -( (1.-integ%alpha_f)*integ%gamma - integ%beta )/integ%beta
     b3 = -( integ%gamma - 2*integ%beta)*(1.-integ%alpha_f)*dt/(2.*integ%beta)
     ! Clear the dyn_rhs container
     body%dyn_rhs(1:neq) = 0.
     ! rhs = -M*vec
     vec(1:neq) = a1*(body%qn(1:neq) - body%qi(1:neq)) + a2*body%qdi(1:neq) + a3*body%qddi(1:neq)
     call DGEMV_sparse( body%neq, body%qq_nnz, body%M, &
                        body%qq_IA, body%qq_JA,        &
                        -1.0, vec, 1.0,                &  
                        body%dyn_rhs )

     ! if there are applied kinematics, then 
     ! $rhs = rhs - M_{qv}*\ddot{v}$
     if( allocated( body%restraints_surf ) ) then

        vec(1:ng) = (1.-integ%alpha_m)*body%qddn(ng1:ndofs) + integ%alpha_m* body%qddi(ng1:ndofs)

        call DGEMV_sparse(body%ng,body%qv_nnz,body%Mqv,  &
                          body%qv_IA,body%qv_JA,         &
                          -1.0, vec, 1.0,                &  
                          body%dyn_rhs )

        if( body%damping_flag == SM_TRUE ) then
           vec(1:ng) = (1.-integ%alpha_f)*body%qdn(ng1:ndofs) + integ%alpha_f*body%qdi(ng1:ndofs)
           call DGEMV_sparse(body%ng,body%qv_nnz,body%Dampqv,  &
                             body%qv_IA,body%qv_JA,            &
                             -1.0, vec, 1.0,                   & 
                             body%dyn_rhs )
        end if
     end if

     ! If there is Prop. Damping
     if( body%damping_flag == SM_TRUE ) then
        vec(1:neq) = b1*(body%qn(1:neq) - body%qi(1:neq)) + b2*body%qdi(1:neq) + b3*body%qddi(1:neq)
        call DGEMV_sparse(body%neq,body%qq_nnz,body%damp, &
             body%qq_IA, body%qq_JA,        &
             -1.0, vec, 1.0,                & 
             body%dyn_rhs )    
     end if

     ! Reduce 
     ! rhs = rhs + Fext - Fint
     body%dyn_rhs(1:neq) = body%dyn_rhs(1:neq) + (1.-integ%alpha_f)*( -body%Qs(1:neq) + body%Hs(1:neq)) &
                           + integ%alpha_f*( -body%Qsn(1:neq) + body%Hsn(1:neq) )

     !-----------------------------------------------------
     ! Build the Total Tangent Stiffness Matrix into K
     !
     nnz = body%qq_nnz
     body%K(1:nnz) = body%K(1:nnz) + a1*body%M(1:nnz)
     if( body%damping_flag == SM_TRUE ) then
        body%K(1:nnz) = body%K(1:nnz) + b1*body%Damp(1:nnz)
     end if

     !-----------------------------------------------------
     ! Solve for the Incremental Displacements dq
     !
     ! Build LU = M, and store in LU_factors

     iopt = 1
     nrhs = 1
     ldb = neq
     call sm_c2fortran_dgssv( iopt, body%neq, body%qq_nnz, nrhs, body%K, &
                               body%qq_ia, body%qq_ja, body%dyn_rhs, ldb, &
                               body%lu_factors, info )
     ! Backsolve into dyn_rhs
     iopt = 2
     call sm_c2fortran_dgssv( iopt, body%neq, body%qq_nnz, nrhs, body%K, &
                              body%qq_ia, body%qq_ja, body%dyn_rhs, ldb, &
                              body%lu_factors, info )
     ! Cleanup
     iopt = 3
     call sm_c2fortran_dgssv( iopt, body%neq, body%qq_nnz, nrhs, body%K, &
                              body%qq_ia, body%qq_ja, body%dyn_rhs, ldb, &
                              body%lu_factors, info )

     ! update qn = qn + dq
     body%qn(1:neq) = body%qn(1:neq) + body%dyn_rhs(1:neq)

     !-----------------------------------------------------
     ! Update qdn and qddn
     b1 = integ%gamma/( integ%beta * dt)
     b2 = -(integ%gamma - integ%beta)/integ%beta
     b3 = -(integ%gamma - 2*integ%beta)/(2.*integ%beta)*dt
     body%qdn(1:neq)  = b1*(body%qn(1:neq) - body%qi(1:neq)) + b2*body%qdi(1:neq) + b3*body%qddi(1:neq)
     b1 = 1./(integ%beta*dt**2)
     b2 = -1./(integ%beta*dt)
     b3 = -(1.-2*integ%beta)/(2.*integ%beta)
     body%qddn(1:neq) = b1*(body%qn(1:neq) - body%qi(1:neq)) + b2*body%qdi(1:neq) + b3*body%qddi(1:neq)


     !-----------------------------------------------------
     ! loop end conditions
     !
     ! compute $\Delta \dot{q}$

     error_qd = maxval(abs(body%dyn_rhs(1:neq))) * integ%gamma/( integ%beta * dt)

#ifdef DEBUG_SOLID
     write(*,'(A,I3,A,G18.10)') '    GA subit ',iter,' error(qd)=',error_qd
#endif

     if( error_qd <= iter_epsilon ) then
        iter_flag = 1

     else if ( iter >= Nmax_iter ) then
        iter_flag = 2
        write(*,'(A,I3,A,I3,A)') '**Problem: body ',ibd,' did not equilibriate in', iter,'iterations'
        
     else
        iter = iter + 1

     end if

  end do

  !-----------------------------------------------------
  ! Convergence Checks for the FSI
  ! Predictor step always never converges
  if( integ%pcflag .eq. SM_PCPREDICTOR) then

     if( integ%predictor_only .eq. SM_TRUE ) then
        ! if there is an override to only have the prediction step.
        ! this is useful if the body's surface is completely prescribed
        integ%pcflag     = SM_PCPREDICTOR
        integ%pcconvflag = SM_PCCONVERGED
        integ%pcerr      = 0.

     else
        ! Carry on to the corrector
        integ%pcflag     = SM_PCCORRECTOR
        integ%pcconvflag = SM_PCNOTCONVERGED
        integ%pcerr      = 1.
     endif

  ! Corrector
  else

     ! Check for convergence
     ! e1 = L1 norm of k+1Y1-kY1 
     error_q  = maxval(abs(body%qn(1:neq) - body%qms(1:neq,1)))
     !error_qd = maxval( )
     e1 = error_q

     ! Load convergence error for the body:        
     integ%pcerr = e1

     if (e1 .le. integ%pcepsilon) then
        integ%pcflag      = SM_PCPREDICTOR
        integ%pcconvflag  = SM_PCCONVERGED
     endif

  end if

  deallocate(vec)

!#ifdef DEBUG_SOLID
!  ! snapshot number #######=9##p##q
!  !     p: step number
!  !     q: subiter number
!  call Driver_getNstep(nstep)
!  iopt = 9000000 + 1000*nstep + sm_integ_subiter
!  call sm_ioWriteSnapshot(ibd,iopt)
!  call sm_ioWriteParticles(ibd,iopt)
!
!#endif


  return

end subroutine sm_GenAlpha_advance

