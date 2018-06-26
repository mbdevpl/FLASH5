!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Contact/FDEMethod/sm_contact_CollisionVolume
!!
!! NAME
!!  gr_sbDistributedForces
!!!
! SYNOPSIS
!!
!!  gr_sbDistributedForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***

subroutine collisionVolume(xyz2,ind,intVol)

#include "constants.h"

	use sm_contact_interface 
	implicit none

!       TO BE RUN IF IND .GE. 3
!       XYZ MUST BE SET TO 0 AFTER EACH TIMESTEP
!       intVol set to zero after each timestep
	
	integer, intent(IN) :: ind
	real, intent(IN),dimension(ind,3):: xyz2
	real, intent(OUT) :: intVol
	
!       local variables
!       Lapack variables
	integer :: INFO,IDUMMY,ILO,IHI,i,mn,j
	real::DUMMY3,DUMMY4,ABNRM,DUMMY1,len_xyz,len_evec1,len_evec2,len_evec3
	real, dimension(3)::eval,WI,SCALE
	real, dimension(10000):: WORK1  ! named it WORK1 because it is called conflicting with constants.h
	real, dimension(3,3) :: evec
	
	real, dimension(1,3)::com,weight1,weight2,weight3
	real, dimension(3,3)::cov
	real, dimension(ind,3):: xyz
	
 
 !
	com (1,:) = (/sum(xyz2(:,1)),sum(xyz2(:,2)),sum(xyz2(:,3))/);

!	do j=1,ind
!	   write(*,*)'point number ', j, xyz2(j,:)
!	enddo

	com=com/dble(ind);      ! Center of mass of collision volume
!		write(*,*)'com ',com(1,1),com(1,2),com(1,3)
 ! Normalize the volume
	do i=1,ind
	   xyz(i,1:3)=xyz2(i,1:3)-com(1,1:3);
	enddo
	
	cov=1./real(ind-1)*matmul(transpose(xyz),xyz); ! Covariance matrix for PCA

	CALL DGEEVX('B','N','V','N',3,cov,3,eval,WI,DUMMY1,3,evec,3,ILO,IHI,&
	&SCALE,ABNRM,DUMMY3,DUMMY4,WORK1,132,IDUMMY,INFO)
	
	if(INFO.ne.0) then
	   write(*,*)'error in calculating eigenvectors'
	endif
	
	
 !       weight matrices store (/ smallest angle between vector of center of mass to intersect node and PCA vector,
 !       index of node for vector with smallest angle between PCA,
 !       magnitude of vector between C.o.M. and node (becomes radius for ellipsoid) /) 
	
	weight1(1,1:3)=(/100.,0.,0./);
	weight2(1,1:3)=(/100.,0.,0./);
	weight3(1,1:3)=(/100.,0.,0./);

	len_evec1=length(evec(1:3,1));
	len_evec2=length(evec(1:3,2));
	len_evec3=length(evec(1:3,3));
	
        do mn=1,ind
	   
!       loop checks each intersection node to see which has the smallest angle between it and each PCA
!       must be 2 checks, as angles slightly larger and smaller must be normalized to zero

	len_xyz=length(xyz(mn,:));
	   
	if(modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,1))/&
	   &((len_xyz)*len_evec1))),PI) .lt. weight1(1,1)) then
	     weight1(1,1)=modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,1))/&
	        &((len_xyz)*len_evec1))),PI);
	     weight1(1,2)=dble(mn);
	elseif(modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,1))/((len_xyz)&
	   &*len_evec1))),PI) .lt. weight1(1,1)) then
	     weight1(1,1)=modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,1))/&
	        &((len_xyz)*len_evec1))),PI);
	     WEIGHT1(1,2)=dble(mn);
	endif
	
	if(modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,2))/&
	   &((len_xyz)*len_evec2))),PI) .lt. weight2(1,1)) then
   	     weight2(1,1)=modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,2))&
	        &/((len_xyz)*len_evec2))),PI);
	     weight2(1,2)=dble(mn);
	elseif(modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,2))/((len_xyz)&
	   &*len_evec2))),PI) .lt. weight2(1,1)) then
	     weight2(1,1)=modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,2))/&
	        &((len_xyz)*len_evec2))),PI);
	     weight2(1,2)=dble(mn);
	endif
	
	if(modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,3))/((len_xyz)&
	   &*len_evec3))),PI) .lt. weight3(1,1)) then
             weight3(1,1)=modulo(PI-abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,3))/&
   	       &((len_xyz)*len_evec3))),PI);
	     weight3(1,2)=dble(mn);
	elseif(modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,3))/((len_xyz)&
	   &*len_evec3))),PI) .lt. weight3(1,1)) then
	     weight3(1,1)=modulo(abs(acos(DOT_PRODUCT(xyz(mn,:),evec(1:3,3))/&
	        &((len_xyz)*len_evec3))),PI);
	     weight3(1,2)=dble(mn);
	endif
	enddo
	
	weight1(1,3)=length(xyz(nint(weight1(1,2)),:));
	weight2(1,3)=length(xyz(nint(weight2(1,2)),:));
	weight3(1,3)=length(xyz(nint(weight3(1,2)),:));		
	
!       Volume of ellipsoid
	intVol=4./3.*PI*(weight1(1,3)*weight2(1,3)*weight3(1,3));
	
	end subroutine collisionVolume
	



