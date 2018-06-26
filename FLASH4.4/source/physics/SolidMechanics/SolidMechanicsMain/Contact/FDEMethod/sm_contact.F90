!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Contact/FDEMethod/sm_contact
!!
!! NAME
!!  sm_Contact
!!
!! SYNOPSIS
!!
!! sm_Contact()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***
#include "SolidMechanics.h"
#include "Flash.h"
  
subroutine  sm_contact(ib,restart_local)

  use SolidMechanics_data, Only : sm_bodyinfo, sm_structure ,sm_nen, &
       sm_NumBodies, sm_Meshme
  use sm_Misc_interface, only: get_Nodal_UVW_qdi,get_Nodal_qn
  use sm_contact_data
  use sm_contact_interface, only : CollisionVolume, detectCollision
  use sm_Verlet_data, only: sm_Verlet_type, sm_Verlet_info
  use driver_data, only: dr_nstep

  implicit none

  interface 
     function length(a)
       real,dimension(3)::a
       real::length
     end function length
  end interface

  ! Argument list
  integer, intent(IN) :: ib
  logical, intent(IN) :: restart_local
  
  ! local variables
  real :: dt_local
  real,dimension(NDIM) :: point,vel
  real :: intersect,intersectOld
  real, save, dimension(MAXNODERBC,NDIM) ::  vel123,XYZ
  real, allocatable :: xyzIntersect(:,:) 
  real :: intVol
  real :: K_s=1000.
  real, dimension(NDIM) :: com1,com2,comm12, F_Coll
  ! integers
  integer :: nen_e,i,j,k,nnp,Nele,kk,ind,m
  integer :: ibd, ibd2,ind1,ind2,jj
  integer, allocatable :: indcies1(:) , indcies2(:) 
  ! data structures
  type(sm_structure),pointer  :: body1,body2 
  type(sm_Verlet_type), pointer :: integ    
  type(contact_structure),pointer :: contact,contact1, contact2
  
 
  ! Get all the specfic data related to body ibd
  body1 => sm_BodyInfo(ib)
  nnp  = body1%nnp;
  contact => sm_contactInfo(ib)
  integ=> sm_Verlet_info(ib)
  dt_local=integ%dt;

  ! Take action only if you are the owner of this Solidbody
  if ( body1 % bodyMaster== sm_MeshMe) then
     select case(body1%bodytype)
        
     case (BODYTYPE_RBC) 
     
     do i=1,NDIM
        body1%comm(i) = sum(body1%qn( (/(j,j=i,NDIM*nnp-(NDIM-i),NDIM)/) ))/real(nnp);
      end do
      body1%comm=body1%comm/real(nnp);
      
      !ind=0;
      do j=1,nnp
         point(:)  = body1%qn((j*3)-2:j*3)    ! xyz of the point on body 1
         vel(:)    = body1%qdms(j*3-2:j*3,1)
         
         intersect=0.
         do i= 1, sm_NumBodies
            if (i.ne.ib) then 
               ind= contact%intersect_count(i);
               body2 => sm_BodyInfo(i)
               do m=1,NDIM

                  body2%comm(m) = sum(body2%qn( (/(jj,jj=m,NDIM*nnp-(NDIM-m),NDIM)/) ))/real(nnp);

#if NDIM==3

#endif
               enddo

               
               Nele  = body2%Nele
               do k=1,Nele            
                  nen_e = sm_nen(body2%eltype(k));         
                  call get_Nodal_qn(k, nen_e, XYZ, body2)  ! positions at beginning of time step
                  call get_Nodal_UVW_qdi(k, nen_e, vel123, body2) ! velocities at the beginning of the timestep
                  ! CHECK IF THIS NODE WAS ALREADY TAGGED IN A PREVIOUS INTERSECTION.
                  !write(*,*) 'before Any Intersect=',intersect,j
                  if (ANY(contact%intersectedNodes(i,:)==j)) intersect=1.
                  !write(*,*) 'After Any Intersect=',intersect,j
                  intersectOld=intersect;
		  
                  
                  
                  if(((point(1) .gt. body1%comm(1)  .and.  point(1) .lt. body2%comm(1)) .or. (point(1) .gt. body2%comm(1)  .and.  point(1) .lt. body1%comm(1))) .or.&
                    &((point(2) .gt. body1%comm(2)  .and.  point(2) .lt. body2%comm(2)) .or. (point(2) .gt. body2%comm(2)  .and.  point(2) .lt. body1%comm(2))) .or.&
                    &((point(3) .gt. body1%comm(3)  .and.  point(3) .lt. body2%comm(3)) .or. (point(3) .gt. body2%comm(3)  .and.  point(3) .lt. body1%comm(3)))) then
                  
                      if(((XYZ(1,1) .gt. body1%comm(1)  .and.  XYZ(1,1) .lt. body2%comm(1)) .or. (XYZ(1,1) .gt. body2%comm(1)  .and.  XYZ(1,1) .lt. body1%comm(1))) .or.&
                        &((XYZ(1,2) .gt. body1%comm(2)  .and.  XYZ(1,2) .lt. body2%comm(2)) .or. (XYZ(1,2) .gt. body2%comm(2)  .and.  XYZ(1,2) .lt. body1%comm(2))) .or.&
                        &((XYZ(1,3) .gt. body1%comm(3)  .and.  XYZ(1,3) .lt. body2%comm(3)) .or. (XYZ(1,3) .gt. body2%comm(3)  .and.  XYZ(1,3) .lt. body1%comm(3)))) then
                          call detectCollision(point,vel,XYZ(1,:),XYZ(2,:),XYZ(3,:), &
                            vel123(1,:),vel123(2,:),vel123(3,:), & 
                            intersect, dt_local) 
!                  if(j==584 .and. dr_nstep .gt. 1) then
!		  	     write(*,*)point(1:3),body1%comm(1),body1%comm(2),body1%comm(3),body2%comm(1),body2%comm(2),body2%comm(3)
!                             stop			     
!	          endif
                          if ((intersect>0.5).and.(intersectOld<0.5)) then
                              ind=ind+1
                              !write(*,*) 'intersect',intersect, ' intersectOld',intersectOld
                              !write(*,*) 'A new node intersects, with total nodes, ',ind, &
                              !  'point ',j, 'on body ',ib
                              contact%intersectedNodes(i,ind)=j;
                              !write(*,*) 'Current nodes from body',ib, 'crossing body,',i,&
                              !contact%intersectedNodes(i,1:ind)
                     
                              exit
                          elseif ((intersect<0.5).and.(intersectOld>0.5)) then  
                              !write(*,*) 'intersect',intersect, ' intersectOld',intersectOld
                              !write(*,*) 'i,',i,'j, ',j,ind
                              !write(*,*) 'A node came out of intersection, new total nodes, ',ind, &
                              ! 'point ',j, 'on body ',ib
                              where (contact%intersectedNodes(i,:)==j)
                                  contact%intersectedNodes(i,:)=contact%intersectedNodes(i,ind)
                              end where
                     
                              contact%intersectedNodes(i,ind)=0
                              !write(*,*) 'Current nodes from body',ib, 'crossing body,',i,&
                              !contact%intersectedNodes(i,1:ind-1)
                              ind=ind-1;
                         end if
                      end if
                  end if
            end do
               


               ! store the number of particles intersected from body ib to body i
               contact%intersect_count(i)=ind;

            end if
         end do
      end do
      
      do ibd = 1, sm_NumBodies-1
         contact1=> sm_ContactInfo(ibd)
         contact1%intersectForce=0.0
         body1 => sm_BodyInfo(ibd) 
         do ibd2=ibd+1,sm_NumBodies
            contact2=> sm_ContactInfo(ibd2)
            contact2%intersectForce=0.0
            body2 => sm_BodyInfo(ibd2) 
            ind1= contact1%intersect_count(ibd2);
            ind2= contact2%intersect_count(ibd) ;
            if (ind1+ind2>3) then
               allocate(xyzIntersect(ind1+ind2,NDIM))
               allocate(indcies1(ind1),indcies2(ind2))
               indcies1=contact1%intersectedNodes(ibd2,1:ind1)
               indcies2=contact2%intersectedNodes(ibd,1:ind2)
               do j=1,ind1+ind2
                  if (j<=ind1) then
                     xyzIntersect(j,:) = body1%qn(indcies1(j)*NDIM-2:indcies1(j)*NDIM)
                  else
                     xyzIntersect(j,:) = body2%qn(indcies2(j-ind1)*NDIM-2:indcies2(j-ind1)*NDIM)
                  end if
               end do
               intVol=0.
               call CollisionVolume(xyzIntersect,ind1+ind2,intVol)
               !write(*,*) 'Calculated Intvol: ',intVol, 'between body',ibd,ibd2
               ! Here comes the COLLISION FORCE
               com1(:) = body1%comm;
               com2(:) = body2%comm
               comm12(:)=(com1-com2)/length(com1-com2); 
               F_Coll=intVol* k_s * Comm12
               
               contact1%intersectForce(ibd2,:) =  F_Coll ; 
               contact2%intersectForce(ibd, :) = -F_Coll ;  
               !write(*,*) 'Collision Force: ', F_Coll

               do j=1,ind1
                  !write(*,*) 'Applying to nodes: ',indcies1(j)*NDIM-2,'to ',indcies1(j)*NDIM
                  body1%Cddms(indcies1(j)*NDIM-2:indcies1(j)*NDIM,1) = F_Coll/real(ind1)
               end do
               do j=1,ind2
                  body2%Cddms(indcies2(j)*NDIM-2:indcies2(j)*NDIM,1) =-F_Coll/real(ind2)
               end do
               deallocate(xyzIntersect,indcies1,indcies2)
               
            end if
         end do
      end do
      
   

   case default

   end select
      
  end if
  return  
   
end subroutine sm_contact
 
