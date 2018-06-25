!  subroutine el02_IntForce_rbc (  )
!
!     
!     Calculates the forces on the nodes particpating in an element
!     
!     Call command from the main subroutine  (p1,p2,p3,area,Areas(i),At,n,c,Vt,ftri)
!


!/////////////////////////////////////////////////////////////////////////////////
!
!
!
      subroutine el02_IntForce_rbc(p1,p2,v1,v2,kp,Lmax_M,sforces,vforces,L)
!
!
!
!     
!////////////////////////////////////////////////////////////////////////////////
      
      use RBCmodule, only :m,ks_M,iteration,debugFlag,modeltype,      &
                           kBT,p_M,visc_coeff,visc,eye,               &
                           gamma_t,gamma_c,Rvisc_Coeff,Coeff1,Coeff2, &
                           dt, viscmodel
      
      implicit none
      real*8  :: p1(3),p2(3),Iij(3),vec(3),Iij_t(3,1)
      real*8  :: sforces(2,3),x,vforces(2,3)
      real*8  :: kp,Lmax_M,L,v
      real(8) :: v1(3),v2(3),vij(3)
      REAL(8) :: FDvisc(1,3),FRvisc(1,3),tr,FRvisc1(3,1)
      REAL(8) :: Wij(3,3),Zij,dWijS(3,3),dWijS_bar(3,3),WWij(9)
      INTEGER :: ii,jj

      vec=p1-p2
      L=sqrt(sum(vec * vec))
      Iij=vec/L
      x=L/Lmax_M
      
      Iij_t=RESHAPE(Iij,(/3,1/));

!     The spring model type ( prefix of SF_ means STRESS FREE )
      
      SELECT CASE (modeltype)
         
      CASE ('SF_FENE_POW')
         
         sforces(1,1:3)=((-ks_M*L/(1-x**2))+(kp/L**real(m)))*Iij
         sforces(2,1:3)=-sforces(1,1:3)
         
      CASE ('SF_WLC_POW')
         
         sforces(1,1:3)=(-(KBT/p_M)*(0.25/((1-x)*(1-x))-0.25+x)+  &
                              (kp/L**real(m)))*Iij
         
         sforces(2,1:3)=-sforces(1,1:3)
      
      CASE ('FENE_POW')
         
         sforces(1,1:3)=((-ks_M*L/(1-x**2))+(kp/L**real(m)))*Iij
         sforces(2,1:3)=-sforces(1,1:3)
         
      CASE ('WLC_POW')
         
         sforces(1,1:3)=(-(KBT/p_M)*((0.25/(1-x**2))-0.25+x)+ &
                              (kp/L**real(m)))*Iij
         
         sforces(2,1:3)=-sforces(1,1:3)   
      
      END SELECT

!C*********************************************************************
!C
!C                              Viscous forces
!C
!C********************************************************************
      
      IF (visc) THEN
         
         vij=v1-v2;
         v=SQRT(vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3))
         
         SELECT CASE (viscmodel)
            
         CASE ('Pep')
            !WRITE(*,*)'Viscous model is  ',viscmodel
            
            IF (v.GT.0) THEN
               
               DO ii=1,3
                  DO jj=ii,3
                     CALL RAND_GEN(Zij)
                     
                     Wij(ii,jj)=Zij-0.5;
                     
                  ENDDO
               ENDDO

!     A matrix of independent Wienner increments 
               
               Wij(2,1)=Wij(1,2);
               Wij(3,1)=Wij(1,3);
               Wij(3,2)=Wij(2,3);

               
               
!     The symmetric random matrix
               
               dWijS=0.5*(Wij+TRANSPOSE(Wij))


!              Trace of the Wij/3

               tr= (Wij(1,1)+ Wij(2,2)+ Wij(3,3))/3.0 ;

!              The traceless symmetric part of the random matrix Wij

               dWijS_bar = dWijS-(tr*Eye);


               FDvisc(1,:)=-(gamma_t*vij)-(gamma_c*dot_product(vij,Iij))*Iij

               FRvisc1= MATMUL((Coeff1*dWijS_bar)+(Coeff2*tr*Eye),Iij_t)


               FRvisc = RESHAPE(FRvisc1,(/1,3/));

               vforces(1,1:3)=FDvisc(1,:)+FRvisc(1,:);
               
               vforces(2,1:3)=-vforces(1,1:3);
            ENDIF
         CASE ('visc1')

            !WRITE(*,*)'Viscous model is  ',viscmodel
            vforces(1,1:3)=-visc_coeff*dot_product(vij,Iij)*Iij;
            
            vforces(2,1:3)=-vforces(1,1:3);
            
         CASE ('visc2')
            
!     WRITE(*,*)'Viscous model is  ',viscmodel
            IF (v > 0) THEN
               
               vforces(1,1:3)=-visc_coeff*vij/v;
               
               vforces(2,1:3)=-vforces(1,1:3);
            ELSE
               vforces=0.0; 
            ENDIF
            
         END SELECT
         
      ELSE
         
         vforces=0.0;
         
         
         
      ENDIF
      
      END SUBROUTINE 
