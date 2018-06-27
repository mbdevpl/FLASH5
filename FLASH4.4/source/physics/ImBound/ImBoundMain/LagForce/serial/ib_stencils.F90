
subroutine ib_stencils(ng,nx,ny,nbd,xb,yb,sb,dsx,dsy, & 
          ielem,jelem,flago,flagi,nxp,nyp,lb,gridflag,&
          np,lpindex,nx1,ny1,nz1,del,coord,bsize)

  use ImBound_data , only : ib_ABODY,ib_nmaxa, ib_stencil, ib_alphax, ib_alphay

  implicit none
#include "Flash.h"
#include "constants.h"

  integer, INTENT(IN) :: ng,nx,ny,nbd,lb,gridflag,nx1,ny1,nz1
  integer, INTENT(INOUT) :: np
  real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa),nxp(ib_nmaxa),nyp(ib_nmaxa)
  real, INTENT(INOUT) :: dsx(ib_nmaxa),dsy(ib_nmaxa)
  integer, INTENT(INOUT) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa)
  integer, INTENT(INOUT) :: flago(nx1,ny1,nz1),flagi(nx1,ny1,nz1),lpindex(ib_nmaxa)
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
   

  ! Local Variables:
  integer :: i,j,k,ii,signflg,count
  integer :: ibd, npoints, indx1, indx2, indy1, indy2
  real    :: cublim(2*NDIM),auxx,auxy,distx,disty,xp,yp
  integer :: icpoint,jcpoint,flagio
     

  integer :: bndx1,bndx2,bndy1,bndy2,bodyflag,ind1,ind2
  integer :: flagw,flage,flags,flagn
  real :: xlower,xupper,ylower,yupper,dx,dy,dxaux,dyaux,xi,yj
  real :: xmin,xmax,ymin,ymax,eps
  real :: aux(ib_nmaxa)


  flago(:,:,:) = 0;
  flagi(:,:,:) = 0;

  dx = del(IAXIS)
  dy = del(JAXIS)

  !write(*,*) nxb,nyb,dx,dy

  eps = 1E-2*MIN(dx,dy)

  bndx1 = ng
  bndx2 = nx     
  bndy1 = ng
  bndy2 = ny     

  if(gridflag .eq. IAXIS) then ! X velocities grid 
     dxaux = 0.0
     dyaux = 0.5*dy

     ! Find block boundaries
     xlower = coord(1) - bsize(1)/2.0 - 1.5*dx
     xupper = coord(1) + bsize(1)/2.0 + 1.5*dx
     ylower = coord(2) - bsize(2)/2.0 - dy
     yupper = coord(2) + bsize(2)/2.0 + dy

  elseif(gridflag .eq. JAXIS) then ! Y velocities grid
     dxaux = 0.5*dx
     dyaux = 0.0
     
     ! Find block boundaries
     xlower = coord(1) - bsize(1)/2.0 - dx
     xupper = coord(1) + bsize(1)/2.0 + dx
     ylower = coord(2) - bsize(2)/2.0 - 1.5*dy
     yupper = coord(2) + bsize(2)/2.0 + 1.5*dy

  endif

  bodyflag = 0
  count = 0

  do ibd=1,nbd

     ! Find bounding box of body
     ind1 = ib_ABODY(ibd) % lb + 1
     ind2 = ib_ABODY(ibd) % lb + ib_ABODY(ibd) % mb

     xmin = MINVAL( xb(ind1:ind2) )
     xmax = MAXVAL( xb(ind1:ind2) )
     ymin = MINVAL( yb(ind1:ind2) )
     ymax = MAXVAL( yb(ind1:ind2) )
   
     ! Check if block is within bounding box of body
     if( xlower .le. xmax .and. xupper .ge. xmin .and. & 
         ylower .le. ymax .and. yupper .ge. ymin ) then

        ! Loop through Lagrangian points:
        do ii = ind1,ind2

           xp = xb(ii);
           yp = yb(ii);

           ! Check if Lagrangian point is within block
           if( xp .ge. xlower .and. xp .le. xupper .and. &
               yp .ge. ylower .and. yp .le. yupper ) then

              count = count + 1
              lpindex(count) = ii

              ! Block on boundary of body
              bodyflag = 1

              ! Find closest point:
              distx = 10.e10;
              disty = 10.e10;
              icpoint = 0
              jcpoint = 0
              do i = bndx1,bndx2
                 xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                      real(i - ng - 1)*dx + dxaux
                 auxx = abs(xp - xi);
                 if (auxx+eps .le. distx) then
                    icpoint = i;
                    distx = auxx;
                 endif
              enddo
    
              do j = bndy1,bndy2
                 yj = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                      real(j - ng - 1)*dy + dyaux
                 auxy = abs(yp - yj);
                 if (auxy+eps .le. disty) then
                    jcpoint = j;
                    disty = auxy;
                 endif
              enddo
    

              ! dsx and dsy:
              dsx(ii) = ib_alphax*dx
              dsy(ii) = ib_alphay*dy

              !write(*,*) 'alphaxy=',ib_alphax,ib_alphay
              !write(*,*) 'dx,dy=',dx,dy,dsx(ii),dsy(ii)

              ! Build stencil structure:
              ! Center North South East West:
              ielem(1:ib_stencil,ii) = (/ icpoint,icpoint,icpoint,     &
                                          icpoint+1,icpoint-1 /)
              jelem(1:ib_stencil,ii) = (/ jcpoint,jcpoint+1,jcpoint-1, &
                                          jcpoint,jcpoint /)

              ! Build flago and flagi matrices:
              if (ib_ABODY(ibd) % memflag .EQ. 1) then ! Solid Body
                 do j = 1, ib_stencil                     
                    xi = coord(IAXIS) - 0.5*bsize(IAXIS) + & 
                         real(ielem(j,ii) - ng - 1)*dx + dxaux
                    yj = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                         real(jelem(j,ii) - ng - 1)*dy + dyaux

                    signflg = int(sign(1.,nxp(ii)*(xi - xp) + &
                                          nyp(ii)*(yj - yp)))              

                    if (signflg .EQ. 1) then
                       flago(ielem(j,ii),jelem(j,ii),1) = -ibd
                    else
                       flago(ielem(j,ii),jelem(j,ii),1) =  ibd
                       flagi(ielem(j,ii),jelem(j,ii),1) = -ibd
                    endif
                 enddo
              end if


           end if
        enddo
     endif
  enddo

  np = count

  if(bodyflag .eq. 1) then

     ! Fill internal points in flago with number ibd:
     do ibd=1,nbd   

        if (ib_ABODY(ibd) % memflag .EQ. 1) then ! Solid Body

           do j = bndy1,bndy2

              do i = bndx1,bndx2
                 if (flago(i-1,j,1) .eq. ibd .and.     &
                     flago(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 elseif (flagi(i-1,j,1) .eq. ibd .and. &
                         flagi(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 endif
              enddo

              do i = bndx2,bndx1,-1
                 if (flago(i+1,j,1) .eq. ibd .and. &
                     flago(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 elseif (flagi(i+1,j,1) .eq. ibd .and. &
                         flagi(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 endif
              enddo

           enddo

           do i = bndx1,bndx2

              do j = bndy1,bndy2
                 if (flago(i,j-1,1) .eq. ibd .and. &
                     flago(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 elseif (flagi(i,j-1,1) .eq. ibd .and. &
                         flagi(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 endif
              enddo

              do j = bndy2,bndy1,-1
                 if (flago(i,j+1,1) .eq. ibd .and. &
                     flago(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 elseif (flagi(i,j+1,1) .eq. ibd .and. &
                         flagi(i,j,1) .eq. 0) then
                    flagi(i,j,1) = ibd
                 endif
              enddo

           enddo

        endif
     enddo

  else

         ! Flag blocks that are completely inside a body
     do ibd=1,nbd

        ! Find bounding box of body
        ind1 = ib_ABODY(ibd) % lb + 1
        ind2 = ib_ABODY(ibd) % lb + ib_ABODY(ibd) % mb

        xmin = MINVAL( xb(ind1:ind2) )
        xmax = MAXVAL( xb(ind1:ind2) )
        ymin = MINVAL( yb(ind1:ind2) )
        ymax = MAXVAL( yb(ind1:ind2) )

        flagw = 0
        flage = 0
        flags = 0
        flagn = 0
            
        if( xlower .gt. xmin .and. xupper .lt. xmax .and. &
            ylower .gt. ymin .and. yupper .lt. ymax ) then

           aux(ind1:ind2) = xlower - xb(ind1:ind2)
               
           if( MAXVAL(aux(ind1:ind2)) .gt. 0.0 ) then

              ! Loop through Eulerian points box:
              do ii = ind1,ind2

                 if(aux(ii) .gt. 0.0) then
                    if( yb(ii) .gt. ylower .and. &
                        yb(ii) .lt. yupper ) then
                       flagw = 1
                       exit
                    endif
                 endif
              enddo
           endif

           if(flagw .ne. 1) goto 123

           aux(ind1:ind2) = xb(ind1:ind2) - xupper

           if( MAXVAL(aux(ind1:ind2)) .gt. 0.0 ) then

              do ii = ind1,ind2

                 if(aux(ii) .gt. 0.0) then
                    if( yb(ii) .gt. ylower .and. &
                        yb(ii) .lt. yupper )  then
                       flage = 1
                       exit
                    endif
                 endif
              enddo
           endif

           if(flage .ne. 1) goto 123

           aux(ind1:ind2) = ylower - yb(ind1:ind2)

           if( MAXVAL(aux(ind1:ind2)) .gt. 0.0 ) then

              do ii = ind1,ind2

                 if(aux(ii) .gt. 0.0) then
                    if( xb(ii) .gt. xlower .and. &
                        xb(ii) .lt. xupper )  then
                       flags = 1
                       exit
                    endif
                 endif
              enddo
           endif

           if(flags .ne. 1) goto 123

           aux(ind1:ind2) = yb(ind1:ind2) - yupper

           if( MAXVAL(aux(ind1:ind2)) .gt. 0.0 ) then

              do ii = ind1,ind2

                 if(aux(ii) .gt. 0.0) then
                    if( xb(ii) .gt. xlower .and. &
                        xb(ii) .lt. xupper )  then
                       flagn = 1
                       exit
                    endif
                 endif
              enddo
           endif

           if(flagn .eq. 1) flagi = ibd
                       
123        continue

        endif
     enddo
  endif

  !write(*,*) nnoda
  !write(*,*) size(xb,DIM=1),size(yb,DIM=1)
  !pause


  !OPEN(UNIT=301,FILE='x.dat',form='formatted')
  !write(301,'(f24.16)') time
  !do i = 1,nx
  !      write(301,'(f24.16)') x(i)
  !enddo
  !close(301)

  !OPEN(UNIT=301,FILE='y.dat',form='formatted')
  !write(301,'(f24.16)') time
  !do i = 1,ny
  !      write(301,'(f24.16)') y(i)
  !enddo
  !close(301)


          
  !OPEN(UNIT=301,FILE='flago.dat',form='formatted')
  !write(301,'(f24.16)') time
  !do i = 1,nx
  !   do j = 1,ny
  !      write(301,'(I8)') flago(i,j,1)
  !   enddo
  !enddo
  !close(301)
  !OPEN(UNIT=301,FILE='flagi.dat',form='formatted')
  !write(301,'(f24.16)') time
  !do i = 1,nx
  !   do j = 1,ny
  !      write(301,'(I8)') flagi(i,j,1)
  !   enddo
  !enddo
  !close(301)

  !pause

  return

End Subroutine ib_stencils
