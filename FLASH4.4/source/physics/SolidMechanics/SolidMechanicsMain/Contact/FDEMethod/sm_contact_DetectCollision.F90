

subroutine detectCollision(point,vel,s1,s2,s3,s_vel1,s_vel2,s_vel3,intersect,dt_local)
  implicit none
  real, dimension(1,3), intent(IN)::point,vel,s1,s2,s3,s_vel1,s_vel2,s_vel3
  real, intent(IN) :: dt_local
  real, intent(INOUT) :: intersect

  real, dimension(1,3) :: norm,a1,a2,a1p,a2p,ndt,point_f
  real, dimension(1,3)::s1_f,s2_f,s3_f,g,pp,s1p,s2p,s3p,d1,d2,vp1,ps1,d1d2,a1d2d1a2,gp
  real, dimension(11)::tv,f1,f2,f12
  
  integer::I,j
  real::dt,b,bdt,xi,zeta,C0,C1,C2,C3,tc,t,tol,error,f,df,tnew,den
 
  interface
     function cross(a,b)
       double precision,dimension(1,3)::a,b,cross
     end function cross
     function length(a)
       double precision,dimension(1,3)::a
       double precision::length
     end function length
  end interface
  
  dt=-dt_local;
  
  point_f=point+vel*dt
  s1_f=s1+s_vel1*dt
  s2_f=s2+s_vel2*dt
  s3_f=s3+s_vel3*dt

  a1=s2_f-s1_f
  a2=s3_f-s1_f
  ndt=cross(a1,a2)
  
  norm=cross(s2-s1,s3-s1)
  
  b=DOT_PRODUCT(norm(1,1:3),(point(1,1:3)-s1(1,1:3)))
  bdt=DOT_PRODUCT(ndt(1,1:3),(point_f(1,1:3)-s1_f(1,1:3)))
  
  
  ! Following is for debugging purposes
  !open(unit=68,file='b_bdt.txt',ACCESS='append')
  !write(68,*)b,bdt,b*bdt,point(1,1:3)-s1(1,1:3),point_f(1,1:3)-s1_f(1,1:3)
  !close(68)
  
  if(b*bdt .le. 0.)then
     
     d1=s_vel2-s_vel1;  ! change of a1 (side vector)
     d2=s_vel3-s_vel1;  ! change of a2
     d1d2=cross(d1,d2);
     vp1=vel-s_vel1;
     ps1=point-s1;

     
     a1d2d1a2=cross(s2-s1,d2)+cross(d1,s3-s1);
     C0=b
     C1=(norm(1,1)*vp1(1,1)+norm(1,2)*vp1(1,2)+norm(1,3)*vp1(1,3))+&
          &DOT_PRODUCT(a1d2d1a2(1,1:3),ps1(1,1:3));
     C2=(d1d2(1,1)*ps1(1,1)+d1d2(1,2)*ps1(1,2)+d1d2(1,3)*ps1(1,3))+&
          &DOT_PRODUCT(a1d2d1a2(1,1:3),vp1(1,1:3));
     C3=(d1d2(1,1)*vp1(1,1)+d1d2(1,2)*vp1(1,2)+d1d2(1,3)*vp1(1,3));
     
     tv=0;
     DO I=1,11
        tv(I)=(I-1)*dt/10.; ! smaller time steps for check
     ENDDO
     
     !writesteps=(/(I,I=firstIter+step,maxIter,step)/)
     f1=C0+C1*tv+C2*tv*tv+C3*tv*tv*tv
     
     f2=(/f1(2:11),f1(11)/);
     
     f12=f1*f2;
     tol=1e-5;
     do j=1,11
        if(f12(j).le.0) then
           tc=tv(j);
           t=f12(j);
        endif
     enddo
     error=100.0;
     
     DO WHILE (error.gt.tol) 
        
        
        f=C0+(C1*t)+(C2*t*t)+(C3*t*t*t);
        df=C1+(2*C2*t)+(3*C3*t*t);
        tnew=t-(f/df);
        
        error=ABS(tnew-t);
        t=tnew;
     ENDDO
     tc=tnew;
     
     pp=point+vel*tc;
     s1p=s1+s_vel1*tc;
     s2p=s2+s_vel2*tc;
     s3p=s3+s_vel3*tc;
     
     gp=pp-s1p;
     
     a1p=s2p-s1p
     a2p=s3p-s1p
     
     den=(length(a1p)*length(a1p))*(length(a2p)*length(a2p))- & 
          DOT_PRODUCT(a1p(1,1:3),a2p(1,1:3))**2;
     
     xi=(DOT_PRODUCT(gp(1,1:3),a1p(1,1:3))*length(a2p)*length(a2p)-DOT_PRODUCT(gp(1,1:3),a2p(1,1:3))&
	  &*DOT_PRODUCT(a1p(1,1:3),a2p(1,1:3)))/den
     zeta=(DOT_PRODUCT(gp(1,1:3),a2p(1,1:3))*length(a1p)*length(a1p)-DOT_PRODUCT(gp(1,1:3),a1p(1,1:3))&
	  &*DOT_PRODUCT(a1p(1,1:3),a2p(1,1:3)))/den
     !if (absintersect( )).and.intersect() 
     if(xi .ge. 0. .and. zeta .ge.  0. .and. xi+zeta .le. 1.) then
        if(intersect .gt. 0.5) then
           intersect=0.;
        else
           intersect=1.;
        endif
     endif
endif

end subroutine detectCollision
