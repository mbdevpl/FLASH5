!     
! File:   Assembly.F
! Author: tim
!
! Created on December 15, 2011, 12:10 PM
! Edited by Hussein Ezzeldin Jan 12 2013.
! hmezz@gwu.edu
! This file is modified to calculate the RBC internal forces assemble the Qs matrix

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_assemble_IntForce_rbc(ibd,flag)

  use SolidMechanics_data, only: sm_nen,  sm_BodyInfo, sm_structure
  use sm_element_interface, only: el05_IntForce_Kirch, el05_IntForce_Biot, & 
                                  el12_IntForce_Kirch, el12_IntForce_Biot, &
                                  el02_IntForce_rbc,el03_IntForce_rbc,     &
                                  el04_IntForce_rbc, getTri_norm
  
  use sm_Misc_interface, only: get_Nodal_UVW_qdm, get_Nodal_qms,           &
                               get_Nodal_UVW_qdi, get_Nodal_qn
  use Driver_data, Only : dr_nstep
  use Driver_interface, only : Driver_abortFlash
  use SolidMechanics_rbc_data, only : fext_M, n_M,stretching_exp
  use Simulation_data, Only: Nplus, Nminus

  implicit none
    
  ! Argument list

  integer, intent(in)    :: ibd ! body number
  integer, intent(in)    :: flag
    
  ! local variables
  type(sm_structure), pointer :: body
  integer,save:: nee, loop1, loop2, idx, nen_e,new_e,ind1,ind2
  integer,save :: e,n1, n2, idx1, modelType,i
  real, save,dimension(MAXNODERBC,NDIM) :: XYZ, Qie, vel,qes
  real, save,dimension (NDIM)  :: norm, cent,norm2,cent2
  real ::  a1,a2
  real :: ks1,ks2,th
  real :: area, aoj, areat,area_ot,volume, vol_ot;
  

  write(*,*) 'This is body no.: ', ibd

  body => sm_BodyInfo(ibd)
  
  ! Initialize  Qs

  body % Qs  = 0.
  
  ! Loop over all the elements to calculate the area and volume terms

  select case (flag)
  case (SM_IOPT_NSTEP)
     
     call getTri_Norm(body%qms(:,1),body%Tri,body%normals,body%centers,body%areas,body%vt,&
          body%nnp,body%nele)
     ! Initialize the forces 
     body%qddms=0.0;
     
  case  (SM_IOPT_NMIDSTEP) 


     call getTri_Norm(body%qn(:)   ,body%Tri,body%normals,body%centers,body%areas,body%vt,&
          body%nnp,body%nele)
     ! Initialize the forces 
    body%qddn=0.0;  
  end select
  
  body%at=sum(body%areas(:))
  write(*,*) 'body%at',body%at

  ! Loop over all the elements
  do e = 1, body%nel

     nen_e = sm_nen(body%eltype(e));  
     nee = NDIM*nen_e;
     new_e=body%new_e(e);

     
     ! get stuff: positions, velocity, forces

     select case (flag)
     case (SM_IOPT_NSTEP)
       
        call get_Nodal_qms(e, nen_e, XYZ, body)  ! positions at beginning of time step
        call get_Nodal_UVW_qdm(e, nen_e, vel, body) ! velocities at the beginning of the timestep

     case  (SM_IOPT_NMIDSTEP) 
        
        call get_Nodal_qn     (e, nen_e, XYZ, body)  ! positions at end of time step
        call get_Nodal_UVW_qdi(e, nen_e, vel, body) ! velocities at end of the timestep
        
     end select
     
     select case( body%eltype(e) )
        
     case( EIGHT_NODE_HEXAHEDRON )
        ! 8-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el05_IntForce_Kirch( qes, XYZ, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el05_IntForce_Biot( qes, XYZ, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        end select
        
     case( TWSEVEN_NODE_HEXAHEDRON )
        ! 27-Node hexahedron
        select case( body%MatType(e) )
        case( MATERIAL_KIRCHHOFF )
           call el12_IntForce_Kirch( qes, XYZ, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        case( MATERIAL_BIOT )
           call el12_IntForce_Biot( qes, XYZ, &
                body%YoungsModulus(e), &
                body%PoissonsRatio(e), Qie)
        end select
        
        
     case (TWO_NODE_LINE)   ! spring forces
        new_e=body%new_e(e);
        modelType=Body%MatType(e)
        select case( body%MatType(e) )
           
        case (MATERIAL_WLC_POW ) 
           
           ks1=body%p_M;
           ks2=body%kp_M(new_e);
           !write(*,*)'Body%Lmax_M(new_e)',Body%Lmax_M(new_e)
           call el02_IntForce_rbc(XYZ(1,:),XYZ(2,:), vel(1,:),vel(2,:),&
                ks1,ks2,Body%Lmax_M(new_e), & 
                modelType, qes)
           

        case ( MATERIAL_WLC_CQ )
           
        case ( MATERIAL_WLC_FENE )
           
        case ( MATERIAL_SF_WLC_POW )
          
           ks1=body%p_M;
           ks2=body%kp_M(new_e);
           

           call el02_IntForce_rbc(XYZ(1,1:NDIM),XYZ(2,1:NDIM), vel(1,1:NDIM),vel(2,1:NDIM),&
                ks1,ks2,Body%Lmax_M(new_e), & 
                modelType, qes)

            if (maxval(abs(qes(1:2,1:3)))>1000) then
               write(*,*) 'ks1, ks2',ks1,ks2
               write(*,*) 'Body%Lmax_M(new_e)',Body%Lmax_M(new_e)
               write(*,*) 'new_e',new_e,body%IEN(1:2,e)
               write(*,*) XYZ(1,:)
               write(*,*) XYZ(2,:)
               write(*,*) qes(1,1:3)
               write(*,*) qes(2,1:3)
               print*, 'Problem in the spring forces'
               stop
           end if
        case ( MATERIAL_SF_WLC_CQ )
           ks1=body%p_M;
           ks2=body%C1;
           
           call el02_IntForce_rbc(XYZ(1,1:NDIM),XYZ(2,1:NDIM), vel(1,1:NDIM),vel(2,1:NDIM),&
                ks1,ks2,Body%Lmax_M(new_e), & 
                modelType, qes)
           
        case (MATERIAL_SF_WLC_FENE)
           
        end select
        
     case(THREE_NODE_TRIANGLE)    ! Area and Volume forces
        new_e=body%new_e(e);        
        norm=body%normals(3*new_e-2:3*new_e);
        cent=body%centers(3*new_e-2:3*new_e);
        
        select case (body%MatType(e))
        
        case (MATERIAL_SF_AREA_VOL)
           !write(*,*) 'Before calling el03'
           area=body%areas(new_e);
           aoj=body%Aoj(new_e);
           areat=body%at;
           area_ot=body%aot;
           volume=body%vt;
           vol_ot=body%vot;
!!$
!!$           call el03_IntForce_rbc(XYZ(1,1:NDIM),XYZ(2,1:NDIM),XYZ(3,1:NDIM), &
!!$                body%areas(new_e),body%Aoj(new_e),body%At,body%Aot,   &
!!$                body%Vt,body%Vot,norm,cent, qes(1:nen_e,1:NDIM))
            call el03_IntForce_rbc(XYZ(1,1:NDIM),XYZ(2,1:NDIM),XYZ(3,1:NDIM), &
                 area,Aoj,Areat,Area_ot,   &
                 Volume,Vol_ot,norm,cent, qes)

           !write(*,*) 'After calling el03'
           if (maxval(abs(qes(1:3,1:3)))>2000) then
              write(*,*)XYZ(1,:)
              write(*,*)XYZ(2,:)
              write(*,*)XYZ(3,:)
              write(*,*)'norm',norm
              write(*,*)'cent',cent
              write(*,*) 'Aoj',body%Aoj(new_e),'Area',body%areas(new_e)
              write(*,*) 'new_e',new_e,body%IEN(1:3,e)
              write(*,*) 'qes'
              write(*,*)qes(1,1:3)
              write(*,*)qes(2,1:3)
              write(*,*)qes(3,1:3)
              print*, 'Problem in the area forces'
              stop
           end if   


        case (MATERIAL_AREA_VOL) ! check this
           call el03_IntForce_rbc(XYZ(1,:),XYZ(2,:),XYZ(3,:), &
                body%areas(e),body%Aoj(e),body%At,body%Aot,   &
                body%Vt,body%Vot,norm,cent, qes)
           
        end select
       
        
     case( FOUR_NODE_QUADRILATERAL )     ! Bending forces

        new_e=body%new_e(e);
        
        !if (new_e==1)print*,"Calculating bending forces"
        select case (body%MatType(e))
        case (MATERIAL_BEND, MATERIAL_SF_BEND)
           
           ind1=body%TRILinks(1,new_e);
           ind2=body%TRILinks(2,new_e);
           !write(*,*)'ind1',ind1,'ind2',ind2
          
           norm =  body%normals(3*ind1-2:3*ind1);
           norm2=  body%normals(3*ind2-2:3*ind2);
           cent =  body%centers(3*ind1-2:3*ind1);
           cent2=  body%centers(3*ind2-2:3*ind2);
           a1   =  body%areas(ind1);
           a2   =  body%areas(ind2);
           
           
           call el04_IntForce_rbc (XYZ(1,:),XYZ(2,:),XYZ(3,:),XYZ(4,:),  &
                norm,norm2,cent,cent2,a1,a2,th,        &
                body%tho, body%cos_tho, body%sin_tho,qes)
           if (maxval(abs(qes(1:4,1:3)))>1000) then
              write(*,*) 'a1,a2',a1,a2 
              write(*,*) 'new_e',new_e,e,ind1,ind2,body%IEN(1:4,e)
              print*, 'Problem in the bending forces'
              stop
           end if
        case default
           call Driver_abortFlash("Error in sm_Assemble_IntForces: unknown FOUR_NODE Material type.")
        end select
        
        
     case default
        call Driver_abortFlash("Error in sm_Assemble_stiff: unknown element type.")
     end select
     
     ! add into global internal force vector
     do loop1 = 1,nen_e
        !write(*,*)'nen_e',nen_e
        idx = body%IEN(loop1,e)
        
        do loop2=1,NDIM
           idx1= body%ID(loop2,idx)
           !write(*,*)'idx',idx,'idx1',idx1
           body%Qs(idx1) = body%Qs(idx1) + qes(loop1,loop2)
           
           if (flag==SM_IOPT_NSTEP) then
              !write(*,*)'1 SM_IOPT_NSTEP' 
              !write(*,*)idx,idx1
             
              body%qddms(idx1,1) = body%qddms(idx1,1) + qes(loop1,loop2)
               !write(*,*)'2 SM_IOPT_NSTEP' 
           elseif (flag==SM_IOPT_NMIDSTEP) then
              !write(*,*) 'SM_IOPT_NMIDSTEP'
              
              body%qddn(idx1)  = body%qddn(idx1) + qes(loop1,loop2)
           
           end if
           
        end do
     end do
     
         
   
  enddo
  !fext_M=50.e-12/n_M;
  if ((stretching_exp==1).and.(flag==SM_IOPT_NSTEP)) then
  ! Positive x 
  body%qddn(Nplus*NDIM-2) = body%qddn(Nplus*NDIM-2)  + fext_M;
  ! Negative x 
  body%qddn(Nminus*NDIM-2)= body%qddn(Nminus*NDIM-2) - fext_M;

  elseif ((flag==SM_IOPT_NMIDSTEP).and.(stretching_exp==1)) then
  ! Positive x 
  body%qddms(Nplus*NDIM-2,1) = body%qddms(Nplus*NDIM-2,1)  + fext_M;
  ! Negative x 
  body%qddms(Nminus*NDIM-2,1)= body%qddms(Nminus*NDIM-2,1) - fext_M;
  end if

end subroutine sm_assemble_IntForce_rbc
