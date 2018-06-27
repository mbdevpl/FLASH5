!  subroutine el13_IntForce_rbc (  )
!
!     
!     Calculates the forces on the nodes particpating in an element
!     
!     Call command from the main subroutine  (p1,p2,p3,area,Areas(i),At,n,c,Vt,ftri)
!



subroutine el13_IntForce_rbc(  )


  !subroutine tri_forces(p1,p2,p3,area,Aoj,At,normal,center,Vt, &
  !                     Aot,Vot,tforces) 

  
  !use RBCmodule, only : NEle,ka,kv,kd,iteration,DebugFlag,NRBC
  
  implicit none
  
  Interface
     function cross_product(a,b)   
       real ::a(3),b(3),cross_product(3)        
     end function cross_product
  end interface
  
  
  real ::beta_a,beta_v,beta_d
      real,dimension(3)::p1,p2,p3
      real,dimension(3)::normal,center
      real,dimension(3)::a32,a13,a21
      real,dimension(9,3)::tforces
      real::alpha
      real,dimension(3)::func1,func2,func3
      real::Aoj,Vt,area,At
      real::Aot,Vot
     

      a32=p3-p2
      a13=p1-p3
      a21=p2-p1
      
      ! functional form
      func1=cross_product(normal,a32)
      func2=cross_product(normal,a13)
      func3=cross_product(normal,a21)


      ! Global area effect
      alpha=-real(ka)*(At-Aot)/(Aot*4*area)
      

      tforces(1,:)=alpha*func1
      tforces(2,:)=alpha*func2
      tforces(3,:)=alpha*func3
       
	! local area constraint forces
      beta_d=-real(kd)*(area-Aoj)/(Aoj*4*area)

      tforces(4,:)=beta_d*func1
      tforces(5,:)=beta_d*func2
      tforces(6,:)=beta_d*func3


!     volume constraint forces
      beta_v=-real(kv)*(Vt-Vot)/(Vot*6.)

      
      
      !write(*,*)'Vt',Vt
      !write(*,*)'Vot',Vot
      
      tforces(7,:)=beta_v*(normal/3.+cross_product(center,a32))
      tforces(8,:)=beta_v*(normal/3.+cross_product(center,a13))
      tforces(9,:)=beta_v*(normal/3.+cross_product(center,a21))
        


end subroutine el13_IntForce_rbc

