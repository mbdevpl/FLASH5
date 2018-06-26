

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine ib_distributedForces(blockID, particleData, vortx, vorty, vortz)

  use Grid_Data, only : gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize

  use ImBound_Data, only : ib_stencil,ib_npol,ib_interp,ib_alphax,ib_alphay,ib_alphaz,ib_dt

  use IncompNS_Data, only : ins_invRe

  use ib_interface, only : ib_stencils,ib_buildABLan

  implicit none
  integer, intent(IN) :: blockID
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortx 
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vorty
  real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortz
  real, intent(INOUT) :: particleData(NPART_PROPS)


  ! Local Variables....
  real :: xp,yp,zp,zL,h,hl,dx,dy,dz,dsx,dsy,dsz,ubdd,vbdd,wbdd,nxp,nyp,nzp
  real, dimension(MDIM) :: xbe,del,coord,bsize

  integer, parameter, dimension(MDIM)      :: grdip = (/ CENTER, CENTER, CENTER /)
  integer, parameter, dimension(MDIM,MDIM) :: grdnu = &
           RESHAPE( (/  FACES, FACES,CENTER, CENTER, FACES, FACES, FACES,CENTER, FACES /), (/MDIM,MDIM /))
                    !    k                   i                  j
  real, parameter, dimension(MDIM)      :: dlip = (/ 0.5, 0.5, 0.5 /)
  real, parameter, dimension(MDIM,MDIM) :: dlnu = &
           RESHAPE( (/ 0., 0., 0.5, 0.5, 0., 0., 0., 0.5, 0. /), (/MDIM,MDIM /))  
                    !   k             i           j
  integer :: presflag,gridfl(MDIM),i,j,idim,gridind,nkij
  integer :: ielem(ib_stencil,MDIM,CONSTANT_TWO)
  real :: dpdn,eps
  real :: delaux(MDIM),xyz_stencil(ib_stencil,MDIM),phile(ib_stencil,NDIM+1)

  real :: zpres,zv(MDIM),nuwx,nuwy,nuwz

  real, pointer, dimension(:,:,:,:) :: solnData

  !integer, parameter :: derivflag = 0 ! Give interpolation functions and their derivatives

  real    :: A(ib_npol,ib_npol), B(ib_npol,ib_stencil)
  real    :: p(ib_npol), gamma(ib_npol), indx(ib_npol)

  real    :: auxx,dist(MDIM),xi,d
  integer :: cpoint(MDIM),icpoint,jcpoint,buildflag
  integer :: bnd(LOW:HIGH,MDIM)
  integer, parameter, dimension(MDIM) :: nf = (/GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC/)

  real :: nu, dt

  nu = ins_invRe
  dt = ib_dt

  ! Get dx,dy
  call Grid_getDeltas(blockID,del)
  call Grid_getBlkCenterCoords(blockID,coord)
  call Grid_getBlkPhysicalSize(blockID,bsize)

  dx=Del(IAXIS)
  dy=Del(JAXIS)
#if NDIM == 3
  dz=Del(KAXIS)
#else
  dz=1.
#endif

  ! Particle data:
  xp  = particleData(POSX_PART_PROP)
  yp  = particleData(POSY_PART_PROP)
  ubdd= particleData(ACCX_PART_PROP)
  vbdd= particleData(ACCY_PART_PROP)
  nxp = particleData(NMLX_PART_PROP)
  nyp = particleData(NMLY_PART_PROP)
  
#if NDIM == 3
  zp  = particleData(POSZ_PART_PROP)
  wbdd= particleData(ACCZ_PART_PROP)
  nzp = particleData(NMLZ_PART_PROP)
#else
  zp  = 0.
  wbdd= 0.
  nzp = 0.
#endif


  ! Function Gimme h:
  ! call ib_normaldistance(dx,dy,dz,h)
  dsx = ib_alphax*dx
  dsy = ib_alphay*dy
#if NDIM == 2
  h   = 0.5*(dsx+dsy)  
  eps = 1.E-2*MIN(dx,dy)  
#elif NDIM == 3
  dsz = ib_alphaz*dz
  h   = 1./3.*(dsx+dsy+dsz)
  eps = 1.E-2*MIN(dx,dy,dz)
#endif

               
  ! External Point Position:
  xbe(IAXIS) = xp + nxp*h
  xbe(JAXIS) = yp + nyp*h
#if NDIM == 3
  xbe(KAXIS) = zp + nzp*h
#else
  xbe(KAXIS) = 0.
#endif
  
  zpres = 0.
  zv(1:MDIM) = 0.
  zL= 0.
  nuwx = 0.
  nuwy = 0.
  nuwz = 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE ! CONSTANT_ZERO = Viscous Forces, CONSTANT_ONE = Pressure

     ! N kij
     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)

     do gridind = 1,nkij

        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdnu(:,gridind)

        ! Auxiliary deltas
        do idim = 1,MDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlnu(idim,gridind))*del(idim)

           bnd(LOW,idim)  = CONSTANT_TWO
           bnd(HIGH,idim) = nf(idim)
           if (gridfl(idim) .eq. CENTER) bnd(HIGH,idim) = bnd(HIGH,idim)- 1 ! Pressure
  
        enddo


        ! Find closest point:
        dist(:)   = 10.e10
        cpoint(:) = 0

        ! Closest point in x,y,z
        do idim = 1,NDIM
           do i = bnd(LOW,idim),bnd(HIGH,idim)
              xi = coord(idim) - 0.5*bsize(idim) + &
                   real(i - NGUARD - 1)*del(idim) + delaux(idim)
              auxx = abs(xbe(idim) - xi);
              if (auxx+eps .le. dist(idim)) then
                 cpoint(idim) = i;
                 dist(idim) = auxx;
              endif
           enddo
        enddo


        icpoint = cpoint(IAXIS)
        jcpoint = cpoint(JAXIS)

        
        ! Obtain Stencil for External Point: Only 2D Linear basis.
        ielem(1:ib_stencil,IAXIS,presflag+1) = (/ icpoint,icpoint,icpoint,     &
                                                  icpoint+1,icpoint-1 /)
        ielem(1:ib_stencil,JAXIS,presflag+1) = (/ jcpoint,jcpoint+1,jcpoint-1, &
                                                  jcpoint,jcpoint /)

        ielem(1:ib_stencil,KAXIS,presflag+1) = 1

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:

        ! Build A and B matrices:
        call ib_buildABLan(ib_stencil,ib_npol,dsx, &
            dsy,xbe(IAXIS),xbe(JAXIS),xyz_stencil(1:ib_stencil,IAXIS), &
            xyz_stencil(1:ib_stencil,JAXIS),ib_interp,A,B,buildflag); 

        ! Obtain gamma coefficients:
        ! Solve for systems coefficients:
        if (ib_interp == 1) then
           p(1) = 1.; p(2) = xp; p(3) = yp;   
        endif

        gamma = p

        call ludcmp(A,ib_npol,ib_npol,indx,d)
        call lubksb(A,ib_npol,ib_npol,indx,gamma)

        ! Obtain Shape functions:
        phile =0.
        do i = 1 , ib_stencil
           do j = 1, ib_npol
              phile(i,1) = phile(i,1) + gamma(j)*B(j,i)
           enddo
        enddo

                          
        if (presflag .eq. CONSTANT_ONE) then         ! Pressure

           ! Point to cell centered Variables:
           call Grid_getBlkPtr(blockID,solnData,CENTER)

           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zpres = zpres + phile(i,1)*solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1));   
           enddo

           ! Release Pointer
           call Grid_releaseBlkPtr(blockID,solnData,CENTER)

           ! Get Pressure approximation at surface marker:
           dpdn = -(ubdd*nxp + vbdd*nyp + wbdd*nzp);
           zL = zpres - dpdn*h;

        else                                         ! Tangent stress
           select case (gridind)
           case(1) 
           ! Only wz:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortz(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwz = nu*zv(gridind)

           case(2)
           ! wx:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortx(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwx = nu*zv(gridind)          

           case(3)
           ! wy:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vorty(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo
           nuwy = nu*zv(gridind)

           end select
           
        end if

     enddo

  enddo

  ! Assign pressure and viscous forces to particleData:
  particleData(PRES_PART_PROP) = zL

  ! Fvisc = -nu (N x W)
  particleData(FXVI_PART_PROP) = (nzp*nuwy-nyp*nuwz)
  particleData(FYVI_PART_PROP) = (nxp*nuwz-nzp*nuwx)
  particleData(FZVI_PART_PROP) = (nyp*nuwx-nxp*nuwy)

  return

end subroutine ib_distributedForces


