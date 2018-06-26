!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_iouttotecplot
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!
!! File written in ASCI format, for speedup of Tecplot load use preplot command. 
!! i.e.: preplot geo.000001.plt geobin.000001.plt
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

#if NDIM==MDIM
#define WRITE_BINARY_GEO
#endif

subroutine sm_iouttotecplot(nstep,time,dt,count)

  use SolidMechanics_Data, only : sm_meshMe,sm_meshComm,sm_NumBodies,sm_BodyInfo
  use sm_assemble_interface, only : sm_assemble_COM
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, only : gr_globalDomain
  use gr_sbData, only : gr_sbBodyInfo
  implicit none

  integer, intent(in) :: nstep,count
  real, intent(in)    :: time, dt

  real :: xi,yi,zi,nmlx,nmly,nmlz,pres,fvx,fvy,fvz,wx,wy,wz
  character(len=6) :: index_ibd,index_name
  integer :: i,ie,ibd

  integer :: ierr

  integer :: x_sgn,y_sgn,z_sgn,nwx,nwy,nwz
  real    :: Lx, Ly, Lz
  logical :: vflag
  real, allocatable, dimension(:) :: vx,vy,vz,vxm,vym,vzm

  integer :: plot_nel

  real, parameter :: eps = 1.e-12

#ifdef WRITE_BINARY_GEO
  character(28) :: filename
  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 visdouble,disdouble
  integer*4 Debug,ZoneType
  character*1 NULLCHR
  integer, allocatable, dimension(:,:) :: icon
#endif



#ifndef WRITE_BINARY_GEO
  !!! Write to Tecplot - Put in IO file:
  write(index_name,"(I6.6)") count
  if (sm_meshMe .eq. MASTER_PE) then
     open(unit=113,file='./IOData/geo.'//index_name, &
          status='replace',form='formatted')
     close(113)
     open(unit=110,file='./IOData/parts.'//index_name,form='formatted',&
          status='replace')
     close(110)
  endif
#else

!-----------------------------------------------------------------------
!                                                         TecPlot set-up
!-----------------------------------------------------------------------
  Debug     = 0
  disdouble = 0
  visdouble = 0
  NULLCHR   = CHAR(0)
!-----------------------------------------------------------------------

  ! Initialize Tecplot Write:
  write(filename,'("./IOData/geo.",i4.4,".",i6.6,".plt")') &
        count,sm_meshMe

#if NDIM == MDIM
  i = TecIni('GEO'//NULLCHR,                            &
             'x y z'//NULLCHR,                          &
             filename//NULLCHR,                         &
             './IOData/'//NULLCHR,                      &
             Debug,visdouble)
#else
  i = TecIni('GEO'//NULLCHR,                            &
             'x y'//NULLCHR,                            &
             filename//NULLCHR,                         &
             './IOData/'//NULLCHR,                      &
             Debug,visdouble)
#endif  /* NDIM */

#endif /* WRITE_BINARY_GEO */

  call MPI_BARRIER(sm_meshComm,ierr)

  !! Loop over bodies and write:
  do ibd = 1,sm_NumBodies
     if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster)then
 

#ifndef WRITE_BINARY_GEO
        ! Write Particles
        open(unit=110,file='./IOData/parts.'//index_name,form='formatted',&
             status='old',position='append')
        write(index_ibd,"(I6.6)") ibd
#if NDIM == MDIM
        write(110,'(A)') 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz", "p", "fvx", "fvy", "fvz", "wx", "wy", "wz"' 
        write(110,*)'ZONE I=',gr_sbBodyInfo(ibd)%totalpart,', T="Parts'//index_ibd//'"',', F=POINT' 
        do i = 1,gr_sbBodyInfo(ibd)%totalpart
           xi  =gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) 
           yi  =gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)
           zi  =gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)
           nmlx=gr_sbBodyInfo(ibd)%particles(NMLX_PART_PROP,i)
           nmly=gr_sbBodyInfo(ibd)%particles(NMLY_PART_PROP,i)
           nmlz=gr_sbBodyInfo(ibd)%particles(NMLZ_PART_PROP,i)
           pres=gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,i)
           fvx =gr_sbBodyInfo(ibd)%particles(FXVI_PART_PROP,i)
           fvy =gr_sbBodyInfo(ibd)%particles(FYVI_PART_PROP,i)
           fvz =gr_sbBodyInfo(ibd)%particles(FZVI_PART_PROP,i)
           wx  =gr_sbBodyInfo(ibd)%particles(VORX_PART_PROP,i)
           wy  =gr_sbBodyInfo(ibd)%particles(VORY_PART_PROP,i)
           wz  =gr_sbBodyInfo(ibd)%particles(VORZ_PART_PROP,i) 
           write(110,'(13g19.10)') xi,yi,zi,nmlx,nmly,nmlz,pres,fvx,fvy,fvz,wx,wy,wz        
        enddo
#else
        ! Write Particles 2D:
        write(110,'(A)') 'VARIABLES = "X", "Y", "nx", "ny", "p", "fvx", "fvy", "wz"'
        write(110,*)'ZONE I=',gr_sbBodyInfo(ibd)%totalpart,', T="Parts'//index_ibd//'"',', F=POINT'
        do i = 1,gr_sbBodyInfo(ibd)%totalpart
           xi  =gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)
           yi  =gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)
           nmlx=gr_sbBodyInfo(ibd)%particles(NMLX_PART_PROP,i)
           nmly=gr_sbBodyInfo(ibd)%particles(NMLY_PART_PROP,i)
           pres=gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,i)
           fvx =gr_sbBodyInfo(ibd)%particles(FXVI_PART_PROP,i)
           fvy =gr_sbBodyInfo(ibd)%particles(FYVI_PART_PROP,i)
           wz  =gr_sbBodyInfo(ibd)%particles(VORZ_PART_PROP,i)
           write(110,'(8g19.10)') xi,yi,nmlx,nmly,pres,fvx,fvy,wz
        enddo
#endif  /* NDIM */
        close(110)      

#else
        ! Dont write Particle Information in Binary mode.
#endif /* WRITE_BINARY_GEO */     

        ! Actual positions:
        allocate(vx(sm_BodyInfo(ibd)%nnp),vy(sm_BodyInfo(ibd)%nnp))
              
        vx(:) =  sm_BodyInfo(ibd)%x(:) + &
                 sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(IAXIS,:))

        vy(:) =  sm_BodyInfo(ibd)%y(:) + &
                 sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,:))

#if NDIM == MDIM
        allocate(vz(sm_BodyInfo(ibd)%nnp))
        vz(:) =  sm_BodyInfo(ibd)%z(:) + &
                 sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,:)) 
#endif /* NDIM */

        write(index_ibd,"(I6.6)") ibd

#ifndef WRITE_BINARY_GEO
        open(unit=113,file='./IOData/geo.'//index_name,form='formatted',&
             status='old',position='append')
#if NDIM == MDIM
        write(113,'(A)') 'VARIABLES = "X", "Y", "Z"'
#else
        write(113,'(A)') 'VARIABLES = "X", "Y"'
#endif  /* NDIM */

#else
        ! No need here.
#endif /* WRITE_BINARY_GEO */



        if (any(sm_BodyInfo(ibd)%OnBoundary(:,:) .eqv. .TRUE.) .and. .true.) then

           Lx = gr_globalDomain(HIGH,IAXIS) - gr_globalDomain(LOW,IAXIS)
           Ly = gr_globalDomain(HIGH,JAXIS) - gr_globalDomain(LOW,JAXIS)
#if NDIM == MDIM
           Lz = gr_globalDomain(HIGH,KAXIS) - gr_globalDomain(LOW,KAXIS)
#endif

           call sm_assemble_COM(ibd)

           ! Shift sign on each direction:
           x_sgn = 0
           if ( sm_BodyInfo(ibd)%OnBoundary(LOW,IAXIS) .and. &
              (.not. sm_BodyInfo(ibd)%OnBoundary(HIGH,IAXIS))) x_sgn = 1

           if ((.not. sm_BodyInfo(ibd)%OnBoundary(LOW,IAXIS)) .and. &
                sm_BodyInfo(ibd)%OnBoundary(HIGH,IAXIS)) x_sgn =-1

           y_sgn = 0
           if ( sm_BodyInfo(ibd)%OnBoundary(LOW,JAXIS) .and. &
              (.not. sm_BodyInfo(ibd)%OnBoundary(HIGH,JAXIS))) y_sgn = 1

           if ((.not. sm_BodyInfo(ibd)%OnBoundary(LOW,JAXIS)) .and. &
                sm_BodyInfo(ibd)%OnBoundary(HIGH,JAXIS)) y_sgn =-1

#if NDIM == MDIM
           z_sgn = 0
           if ( sm_BodyInfo(ibd)%OnBoundary(LOW,KAXIS) .and. &
              (.not. sm_BodyInfo(ibd)%OnBoundary(HIGH,KAXIS))) z_sgn = 1

           if ((.not. sm_BodyInfo(ibd)%OnBoundary(LOW,KAXIS)) .and. &
                sm_BodyInfo(ibd)%OnBoundary(HIGH,KAXIS)) z_sgn =-1
#endif
           ! Number of intervals in each direction:
           nwx = 0
           if (x_sgn .eq. 1) then 
              nwx = int(abs(maxval(vx(:)+eps)-gr_globalDomain(HIGH,IAXIS))/Lx) 
           elseif(x_sgn .eq. -1) then 
              nwx = int((minval(vx(:)-eps)-gr_globalDomain(LOW,IAXIS))/Lx)
           endif

           nwy = 0
           if (y_sgn .eq. 1) then 
              nwy = int(abs(maxval(vy(:)+eps)-gr_globalDomain(HIGH,JAXIS))/Ly) 
           elseif(y_sgn .eq. -1) then 
              nwy = int((minval(vy(:)-eps)-gr_globalDomain(LOW,JAXIS))/Ly)
           endif

#if NDIM == MDIM           
           nwz = 0
           if (z_sgn .eq. 1) then 
              nwz = int(abs(maxval(vz(:)+eps)-gr_globalDomain(HIGH,KAXIS))/Lz) 
           elseif(z_sgn .eq. -1) then 
              nwz = int((minval(vz(:)-eps)-gr_globalDomain(LOW,KAXIS))/Lz) 
           endif
#endif           
           plot_nel = 0
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( TWO_NODE_LINE)

                 plot_nel = plot_nel + 1

              case( THREE_NODE_TRIANGLE )

                 plot_nel = plot_nel + 1
                 
              case( FOUR_NODE_QUADRILATERAL )
                 
                 plot_nel = plot_nel + 1
                 
              case( NINE_NODE_QUADRILATERAL )
                 
                 plot_nel = plot_nel + 4 
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
           enddo

           ! First Body
           ! Here we distinguish among triangle or quadrilateral elements:  

#ifndef  WRITE_BINARY_GEO

#if NDIM == MDIM
           write(113,*)                                                &
           'ZONE T=Body'//index_ibd//', N=',                           &
           sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
           ', F = FEPOINT, ET = QUADRILATERAL'
           write(113,*)'DT = (SINGLE SINGLE SINGLE)'
           do i = 1,sm_BodyInfo(ibd)%nnp
              xi= vx(i) + real(x_sgn*nwx)*Lx
              yi= vy(i) + real(y_sgn*nwy)*Ly
              zi= vz(i) + real(z_sgn*nwz)*Lz
              write(113,'(3F16.11)') xi,yi,zi
           enddo
#else
           write(113,*)                                                &
           'ZONE T=Body'//index_ibd//', N=',                           &
           sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
           ', F = FEPOINT, ET = QUADRILATERAL'
           write(113,*)'DT = (SINGLE SINGLE)'
           do i = 1,sm_BodyInfo(ibd)%nnp
              xi= vx(i) + real(x_sgn*nwx)*Lx
              yi= vy(i) + real(y_sgn*nwy)*Ly
              write(113,'(2F16.11)') xi,yi
           enddo
#endif /* NDIM */


           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))
                 
              case( TWO_NODE_LINE ) ! Used here in context of rigid bodies

                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%Borigin_node,       &
                      sm_BodyInfo(ibd)%Borigin_node


              case( THREE_NODE_TRIANGLE )
                 
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i)
                 
              case( FOUR_NODE_QUADRILATERAL )
                 
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(4,i)
                 
              case( NINE_NODE_QUADRILATERAL )
                 
                 ! Break in four linear quads:
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(8,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(8,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(7,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(4,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(7,i)

              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
              
           enddo

#else

           ! Here TecZone
           i = TecZne(                                       &
                   'ZONE T=Body'//index_ibd//NULLCHR,        &
                    sm_BodyInfo(ibd)%nnp,                    & 
                    plot_nel,                                &
                    1,                                       & ! Quad
                    'FEBLOCK'//NULLCHR,NULLCHR)


           allocate(vxm(sm_BodyInfo(ibd)%nnp))
           allocate(vym(sm_BodyInfo(ibd)%nnp))
           allocate(vzm(sm_BodyInfo(ibd)%nnp))

           ! Dump node positions:
           vxm(:)= vx(:) + real(x_sgn*nwx)*Lx
           vym(:)= vy(:) + real(y_sgn*nwy)*Ly
           vzm(:)= vz(:) + real(z_sgn*nwz)*Lz
           i = TecDat(sm_BodyInfo(ibd)%nnp,vxm(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vym(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vzm(1),1)

           allocate(icon(4,plot_nel))

           ie = 0;
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( THREE_NODE_TRIANGLE )

                 ie = ie + 1   
                 icon(1:3,ie) = sm_BodyInfo(ibd)%ws_IEN(1:3,i)
                 icon(4,ie)   = icon(3,ie)
                 
              case( FOUR_NODE_QUADRILATERAL )


                 ie = ie + 1   
                 icon(1:4,ie) = sm_BodyInfo(ibd)%ws_IEN(1:4,i)
                 
              case( NINE_NODE_QUADRILATERAL )

                 ! Break in four linear quads:
                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(1,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(2,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(4,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(3,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
            
           enddo

           ! Dump Connectivities
           i = TecNod(icon)

           deallocate(icon)
           deallocate(vxm,vym,vzm)

#endif /* WRITE_BINARY_GEO */ 

           ! Test to write mirror body, only if body is crossing a boundary:
           vflag = .false.

           ! Global Domain low
           if( any( ((vx(:)+real(x_sgn*nwx)*Lx-gr_globalDomain(LOW,IAXIS)) .le. 0.) ).and. &
               any( ((vx(:)+real(x_sgn*nwx)*Lx-gr_globalDomain(LOW,IAXIS)) .ge. 0.) )) vflag = .true.

           if( any( ((vy(:)+real(y_sgn*nwy)*Ly-gr_globalDomain(LOW,JAXIS)) .le. 0.) ).and. &
               any( ((vy(:)+real(y_sgn*nwy)*Ly-gr_globalDomain(LOW,JAXIS)) .ge. 0.) )) vflag = .true.

#if NDIM == MDIM           
           if( any( ((vz(:)+real(z_sgn*nwz)*Lz-gr_globalDomain(LOW,KAXIS)) .le. 0.) ).and. &
               any( ((vz(:)+real(z_sgn*nwz)*Lz-gr_globalDomain(LOW,KAXIS)) .ge. 0.) )) vflag = .true. 
#endif /* NDIM */           

           ! Global Domain high
           if( any( ((vx(:)+real(x_sgn*nwx)*Lx-gr_globalDomain(HIGH,IAXIS)) .le. 0.) ).and. &
               any( ((vx(:)+real(x_sgn*nwx)*Lx-gr_globalDomain(HIGH,IAXIS)) .ge. 0.) )) vflag = .true.

           if( any( ((vy(:)+real(y_sgn*nwy)*Ly-gr_globalDomain(HIGH,JAXIS)) .le. 0.) ).and. &
               any( ((vy(:)+real(y_sgn*nwy)*Ly-gr_globalDomain(HIGH,JAXIS)) .ge. 0.) )) vflag = .true.

#if NDIM == MDIM 
           if( any( ((vz(:)+real(z_sgn*nwz)*Lz-gr_globalDomain(HIGH,KAXIS)) .le. 0.) ).and. &
               any( ((vz(:)+real(z_sgn*nwz)*Lz-gr_globalDomain(HIGH,KAXIS)) .ge. 0.) )) vflag = .true. 
#endif /* NDIM */

           if(vflag .and. .true.)then

              plot_nel = 0
              do i=1,sm_BodyInfo(ibd)%ws_nel

                 select case(sm_bodyInfo(ibd)%ws_eltype(i))

                 case( TWO_NODE_LINE )
 
                    plot_nel = plot_nel + 1

                 case( THREE_NODE_TRIANGLE )

                    plot_nel = plot_nel + 1
                 
                 case( FOUR_NODE_QUADRILATERAL )
                 
                    plot_nel = plot_nel + 1
                 
                 case( NINE_NODE_QUADRILATERAL )
                 
                    plot_nel = plot_nel + 4 
                 
                 case default
                    call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
                 end select
              enddo

              ! Mirror Body
#ifndef  WRITE_BINARY_GEO

#if NDIM == MDIM
              write(113,*)                                                &
              'ZONE T=MirrorBody'//index_ibd//', N=',                     &
              sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
              ', F = FEPOINT, ET = QUADRILATERAL'
              write(113,*)'DT = (SINGLE SINGLE SINGLE)'
              do i = 1,sm_BodyInfo(ibd)%nnp
                 xi= vx(i) + real(x_sgn*(nwx+1))*Lx
                 yi= vy(i) + real(y_sgn*(nwy+1))*Ly
                 zi= vz(i) + real(z_sgn*(nwz+1))*Lz
                 write(113,'(3F16.11)') xi,yi,zi
              enddo
#else
              write(113,*)                                                &
              'ZONE T=MirrorBody'//index_ibd//', N=',                     &
              sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
              ', F = FEPOINT, ET = QUADRILATERAL'
              write(113,*)'DT = (SINGLE SINGLE)'
              do i = 1,sm_BodyInfo(ibd)%nnp
                 xi= vx(i) + real(x_sgn*(nwx+1))*Lx
                 yi= vy(i) + real(y_sgn*(nwy+1))*Ly
                 write(113,'(2F16.11)') xi,yi
              enddo
#endif /* NDIM */

              do i=1,sm_BodyInfo(ibd)%ws_nel

                 select case(sm_bodyInfo(ibd)%ws_eltype(i))
  
                 case( TWO_NODE_LINE ) ! Used here in context of rigid bodies
  
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                         sm_BodyInfo(ibd)%Borigin_node,       &
                         sm_BodyInfo(ibd)%Borigin_node


                 case( THREE_NODE_TRIANGLE )

                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(3,i)

                 case( FOUR_NODE_QUADRILATERAL )
                    
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(4,i)

                 case( NINE_NODE_QUADRILATERAL )
                    
                    ! Break in four linear quads:
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(8,i)
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(9,i)
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(8,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(7,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(4,i)
                    write(113,'(4(I8))')                      &
                         sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                         sm_BodyInfo(ibd)%ws_IEN(7,i)

                 case default
                    call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
                 end select

              enddo

#else

           ! Here TecZone
           i = TecZne(                                       &
                   'ZONE T=MirrorBody'//index_ibd//NULLCHR,        &
                    sm_BodyInfo(ibd)%nnp,                    & 
                    plot_nel,                                &
                    1,                                       & ! Quad
                    'FEBLOCK'//NULLCHR,NULLCHR)


           allocate(vxm(sm_BodyInfo(ibd)%nnp))
           allocate(vym(sm_BodyInfo(ibd)%nnp))
           allocate(vzm(sm_BodyInfo(ibd)%nnp))

           ! Dump node positions:
           vxm(:)= vx(:) + real(x_sgn*(nwx+1))*Lx
           vym(:)= vy(:) + real(y_sgn*(nwy+1))*Ly
           vzm(:)= vz(:) + real(z_sgn*(nwz+1))*Lz
           i = TecDat(sm_BodyInfo(ibd)%nnp,vxm(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vym(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vzm(1),1)

           allocate(icon(4,plot_nel))

           ie = 0;
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( THREE_NODE_TRIANGLE )

                 ie = ie + 1   
                 icon(1:3,ie) = sm_BodyInfo(ibd)%ws_IEN(1:3,i)
                 icon(4,ie)   = icon(3,ie)
                 
              case( FOUR_NODE_QUADRILATERAL )


                 ie = ie + 1   
                 icon(1:4,ie) = sm_BodyInfo(ibd)%ws_IEN(1:4,i)
                 
              case( NINE_NODE_QUADRILATERAL )

                 ! Break in four linear quads:
                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(1,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(2,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(4,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(3,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
            
           enddo

           ! Dump Connectivities
           i = TecNod(icon)

           deallocate(icon)
           deallocate(vxm,vym,vzm)

#endif /* WRITE_BINARY_GEO */

           endif !vflag



        else  ! BODY NOT ON BOUNDARY


           plot_nel = 0
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( TWO_NODE_LINE )

                 plot_nel = plot_nel + 1

              case( THREE_NODE_TRIANGLE )

                 plot_nel = plot_nel + 1
                 
              case( FOUR_NODE_QUADRILATERAL )
                 
                 plot_nel = plot_nel + 1
                 
              case( NINE_NODE_QUADRILATERAL )
                 
                 plot_nel = plot_nel + 4 
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
           enddo

#ifndef  WRITE_BINARY_GEO

#if NDIM == MDIM
           write(113,*)                                                     &
                'ZONE T=Body'//index_ibd//', N=',                           &
                sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
                ', F = FEPOINT, ET = QUADRILATERAL'
           write(113,*)'DT = (SINGLE SINGLE SINGLE)'
           do i = 1,sm_BodyInfo(ibd)%nnp
              xi= vx(i)
              yi= vy(i)
              zi= vz(i)
              write(113,'(3F16.11)') xi,yi,zi
           enddo
#else
           write(113,*)                                                     &
                'ZONE T=Body'//index_ibd//', N=',                           &
                sm_BodyInfo(ibd)%nnp,', E=',plot_nel,                       &
                ', F = FEPOINT, ET = QUADRILATERAL'
           write(113,*)'DT = (SINGLE SINGLE)'
           do i = 1,sm_BodyInfo(ibd)%nnp
              xi= vx(i)
              yi= vy(i)
              write(113,'(2F16.11)') xi,yi
           enddo
#endif
 
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( TWO_NODE_LINE ) ! Used here in context of rigid bodies

                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%Borigin_node,       &
                      sm_BodyInfo(ibd)%Borigin_node

              case( THREE_NODE_TRIANGLE )
   
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i)
                 
              case( FOUR_NODE_QUADRILATERAL )

                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(4,i)
                 
              case( NINE_NODE_QUADRILATERAL )

                 ! Break in four linear quads:
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(8,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(5,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(8,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(7,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(4,i)
                 write(113,'(4(I8))')                      &
                      sm_BodyInfo(ibd)%ws_IEN(9,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(6,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(3,i),        &
                      sm_BodyInfo(ibd)%ws_IEN(7,i)
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
            
           enddo
#else      
           ! Here TecZone
           i = TecZne(                                       &
                   'ZONE T=Body'//index_ibd//NULLCHR,        &
                    sm_BodyInfo(ibd)%nnp,                    & 
                    plot_nel,                                &
                    1,                                       & ! Quad
                    'FEBLOCK'//NULLCHR,NULLCHR)


           ! Dump node positions:
           i = TecDat(sm_BodyInfo(ibd)%nnp,vx(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vy(1),1)
           i = TecDat(sm_BodyInfo(ibd)%nnp,vz(1),1)

           allocate(icon(4,plot_nel))

           ie = 0;
           do i=1,sm_BodyInfo(ibd)%ws_nel

              select case(sm_bodyInfo(ibd)%ws_eltype(i))

              case( THREE_NODE_TRIANGLE )

                 ie = ie + 1   
                 icon(1:3,ie) = sm_BodyInfo(ibd)%ws_IEN(1:3,i)
                 icon(4,ie)   = icon(3,ie)
                 
              case( FOUR_NODE_QUADRILATERAL )


                 ie = ie + 1   
                 icon(1:4,ie) = sm_BodyInfo(ibd)%ws_IEN(1:4,i)
                 
              case( NINE_NODE_QUADRILATERAL )

                 ! Break in four linear quads:
                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(1,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(5,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(2,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(8,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(4,i)

                 ie = ie + 1
                 icon(1,ie) = sm_BodyInfo(ibd)%ws_IEN(9,i)
                 icon(2,ie) = sm_BodyInfo(ibd)%ws_IEN(6,i)
                 icon(3,ie) = sm_BodyInfo(ibd)%ws_IEN(3,i)
                 icon(4,ie) = sm_BodyInfo(ibd)%ws_IEN(7,i)
                 
              case default
                 call Driver_abortFlash("sm_outtotecplot: Element type not recognized.")
              end select
            
           enddo

           ! Dump Connectivities
           i = TecNod(icon)

           deallocate(icon)

#endif /* WRITE_BINARY_GEO */

        endif

#ifndef  WRITE_BINARY_GEO
        close(113)
#endif

        deallocate(vx,vy)
#if NDIM == MDIM
        deallocate(vz)
#endif /* NDIM */



     endif
#ifndef  WRITE_BINARY_GEO
     call MPI_BARRIER(sm_meshComm,ierr)
#endif
  enddo

#ifndef  WRITE_BINARY_GEO
  !!!! -----------
  ! Include time-nstep stamp:
  open(unit=113,file='./IOData/geo.'//index_name,form='formatted',&
       status='old',position='append')
  write(113,'(A,I6,A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T="n=',nstep,' T=',time,'"' 
  close(113)
#else
  i = TecEnd()
#endif

  return

end subroutine sm_iouttotecplot
