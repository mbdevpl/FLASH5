
!  subroutine el02_IntForce_rbc (  )
!
!     
!     Calculates the forces on the nodes particpating in an element
!
!
!
subroutine el02_IntForce_rbc(p1,p2, &
     v1,v2, ks1,ks2,Lmax_M,modelType, &
     sforces)
       
#include "Flash.h"
#include "SolidMechanics.h"
  
  use SolidMechanics_rbc_data, only : sm_rbc_Eye, sm_rbc_gamma_t,sm_rbc_gamma_c,  &
       dpd_c1,dpd_c2, sm_rbc_viscmodel, KBT, m_exp, sm_rbc_viscCoeff
  
  
  implicit none
  real,intent(IN)  :: p1(NDIM),p2(NDIM)
  real,dimension(MAXNODERBC,NDIM),intent(OUT) :: sforces
  real,intent(IN)  :: ks1,ks2,Lmax_M,v1(NDIM),v2(NDIM)
  integer, intent(IN) :: modelType
  
  
  ! local arguments
  real :: Iij(NDIM),vec(NDIM),Iij_t(NDIM,CONSTANT_ONE)
  real :: x, L     
  real :: v,vij(NDIM)
  REAL :: FDvisc(1,3),FRvisc(1,3),tr,FRvisc1(3,1)
  REAL :: Wij(3,3),Zij,dWijS(3,3),dWijS_bar(3,3),WWij(9)
  integer :: ii,jj
  real,dimension(MAXNODERBC,NDIM) :: vforces
  
  ! Calculate the distance between the particles
  vec=p1-p2;
  L=sqrt(sum(vec * vec))
  
  Iij=vec/L;
  x=L/Lmax_M;
  !write(*,*) 'L',L,'Lmax_M',Lmax_M
  
  
  !     The spring model type ( prefix of SF_ means STRESS FREE )
  
  SELECT CASE (modeltype)
     
  CASE (MATERIAL_SF_FENE_POW)
     
     
     sforces(1,1:3)=((-ks1*L/(1-x**2))+(ks2/L**real(m_exp)))*Iij
     sforces(2,1:3)=-sforces(1,1:3)
     
  CASE (MATERIAL_SF_WLC_POW)
     
     sforces(1,1:3)=(-(KBT/ks1)*(0.25/((1.-x)*(1.-x))-0.25+x)+(ks2/L**real(m_exp)))*Iij
     sforces(2,1:3)=-1.*sforces(1,1:3)
     if (maxval(abs(sforces)) >1000) then
        !write(*,*) 'ks1',ks1,'x',x,'ks2',ks2,'L',L,'Iij',Iij
        !write(*,*) 'problem with spring forces'
     !stop
  end if
  CASE (MATERIAL_FENE_POW)
     
     sforces(1,1:3)=((-ks1*L/(1-x**2))+(ks2/L**real(m_exp)))*Iij
     sforces(2,1:3)=-sforces(1,1:3)
     
  CASE (MATERIAL_WLC_POW)
     
     sforces(1,1:3)=(-(KBT/ks1)*((0.25/(1-x**2))-0.25+x)+(ks2/L**real(m_exp)))*Iij
     sforces(2,1:3)=-sforces(1,1:3)   
     
  END SELECT
  
  !C*********************************************************************
  !C
  !C                              Viscous forces
  !C
  !C********************************************************************
  
  ! Viscous forces initialized in case of Non_viscous
  
  vforces(:,:)=0.0;
  
  Iij_t=RESHAPE(Iij,(/3,1/));
  vij=v1-v2;
  v=sqrt(vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3))
  
  
  if (v.gt.0) then
     !write(*,*)'Viscous model is ',sm_rbc_viscmodel
     select case (sm_rbc_viscmodel)
        
     case (Pep_viscModel)
        
        
        !   IF (v.gt.0) THEN
        
        do ii=1,3
           DO jj=ii,3
              CALL RAND_GEN(Zij)
              
              Wij(ii,jj)=Zij-0.5;
              
           END DO
        END DO
        
        !     A matrix of independent Wienner increments 
        
        Wij(2,1)=Wij(1,2);
        Wij(3,1)=Wij(1,3);
        Wij(3,2)=Wij(2,3);
        
        !     The symmetric random matrix
        
        dWijS=0.5*(Wij+TRANSPOSE(Wij))
        !              Trace of the Wij/3
        tr= (Wij(1,1)+ Wij(2,2)+ Wij(3,3))/3.0 ;
        
        !              The traceless symmetric part of the random matrix Wij
        
        dWijS_bar = dWijS-(tr*sm_rbc_Eye);
        FDvisc(1,:)=-(sm_rbc_gamma_t*vij)-(sm_rbc_gamma_c*dot_product(vij,Iij))*Iij
        FRvisc1= MATMUL((dpd_c1*dWijS_bar)+(dpd_c2*tr*sm_rbc_Eye),Iij_t);
        FRvisc = RESHAPE(FRvisc1,(/1,3/));
        
        vforces(1,1:3)=FDvisc(1,:)+FRvisc(1,:);
        vforces(2,1:3)=-vforces(1,1:3);
        
     case (visc_relVel1)
        !write(*,*)'sm_rbc_viscCoeff',sm_rbc_viscCoeff
        vforces(1,1:3)=-sm_rbc_viscCoeff*dot_product(vij,Iij)*Iij;
        vforces(2,1:3)=-vforces(1,1:3);
        
     case (visc_relVel2)
        !write(*,*)'sm_rbc_viscCoeff',sm_rbc_viscCoeff
        vforces(1,1:3)=-sm_rbc_viscCoeff*vij/v;
        vforces(2,1:3)=-vforces(1,1:3);
        
        
     end select
     
  end if

  ! Add the spring and viscous forces 
  sforces(1:2,1:NDIM)=sforces(1:2,1:NDIM)+vforces(1:2,1:NDIM);
  !write(*,*) 'Forces magnitude=    ',sqrt(sum(sforces(1,1:3)*Sforces(1,1:3))) 
  return;

end subroutine el02_IntForce_rbc


