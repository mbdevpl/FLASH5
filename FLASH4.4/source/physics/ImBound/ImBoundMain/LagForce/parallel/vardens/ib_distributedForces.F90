

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine ib_distributedForces(blockID, particleData, vortx, vorty, vortz)

  use Grid_Data, only : gr_meshMe

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &
                             Grid_getDeltas, Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize

  use ImBound_Data, only : ib_stencil,ib_alphax,ib_alphay,ib_alphaz,ib_dt

  use ib_interface, only : ib_stencils,ib_getInterpFunc

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
  integer :: presflag,gridfl(MDIM),i,idim,gridind,nkij
  integer :: ielem(ib_stencil,MDIM,CONSTANT_TWO)
  real :: dpdn,eps
  real :: delaux(MDIM),xyz_stencil(ib_stencil,MDIM),phile(ib_stencil,NDIM+1)

  real :: zpres,zv(MDIM),muwx,muwy,muwz
  integer :: i_ind,j_ind,k_ind

  real, pointer, dimension(:,:,:,:) :: solnData

  integer, parameter :: derivflag = 0 ! Give interpolation functions and their derivatives

  real ::mu, dt

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
  muwx = 0.
  muwy = 0.
  muwz = 0.
  do presflag = CONSTANT_ZERO,CONSTANT_ONE

     ! N kij
     nkij = 1 + (1-presflag)*(2*NDIM-MDIM-1)

     do gridind = 1,nkij

        ! Define Grids in case of vorticity and pressure:
        gridfl(:) = presflag*grdip(:) + (1-presflag)*grdnu(:,gridind)

        ! Auxiliary deltas
        do idim = 1,MDIM
           delaux(idim) = (real(presflag)*dlip(idim)+real(1-presflag)*dlnu(idim,gridind))*del(idim)
        enddo

        ! Obtain Stencil for External Point:
        call ib_stencils(xbe,gridfl,del,coord,bsize,   & 
                         ielem(:,:,presflag+1),hl,COMPUTE_FORCES)

        ! Compute shape functions
        ! Positions of points on the stencil:
        xyz_stencil(1:ib_stencil,1:MDIM) = 0. 
        do idim = 1,NDIM
           xyz_stencil(1:ib_stencil,idim) = coord(idim) - 0.5*bsize(idim) + &
                real(ielem(1:ib_stencil,idim,presflag+1) - NGUARD - 1)*del(idim) + delaux(idim) 
        enddo

        ! Get interpolation functions:
        call ib_getInterpFunc(xbe,xyz_stencil,del,derivflag,phile)

        ! Point to cell centered Variables:
        call Grid_getBlkPtr(blockID,solnData,CENTER)

        if (presflag .eq. CONSTANT_ONE) then         ! Pressure


           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zpres = zpres + phile(i,1)*solnData(PRES_VAR,ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1));   
           enddo

           ! Get Pressure approximation at surface marker:
           dpdn = -(ubdd*nxp + vbdd*nyp + wbdd*nzp);
           zL = zpres - dpdn*h;

        else                                         ! Tangent stress

           ! Closest point location
           i = 1
           i_ind = ielem(i,IAXIS,presflag+1)
           j_ind = ielem(i,JAXIS,presflag+1)           
           k_ind = ielem(i,KAXIS,presflag+1)

           select case (gridind)
           case(1)            
           ! Only wz:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortz(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo

           mu = 0.25*(solnData(VISC_VAR,i_ind-1,j_ind-1,k_ind) + &
                      solnData(VISC_VAR,i_ind-1,j_ind,k_ind)   + &
                      solnData(VISC_VAR,i_ind,j_ind,k_ind)     + &
                      solnData(VISC_VAR,i_ind,j_ind-1,k_ind))

           muwz = mu*zv(gridind)

           case(2)
           ! wx:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vortx(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo

           mu = 0.25*(solnData(VISC_VAR,i_ind,j_ind-1,k_ind-1) + &
                      solnData(VISC_VAR,i_ind,j_ind-1,k_ind)   + &
                      solnData(VISC_VAR,i_ind,j_ind,k_ind)     + &
                      solnData(VISC_VAR,i_ind,j_ind,k_ind-1))

           muwx = mu*zv(gridind)          

           case(3)
           ! wy:
           ! Value of the function in xbe:
           do i = 1 , ib_stencil      
              zv(gridind) = zv(gridind) + phile(i,1)*vorty(ielem(i,IAXIS,presflag+1), &
                                                           ielem(i,JAXIS,presflag+1), &
                                                           ielem(i,KAXIS,presflag+1)); 
           enddo

           mu = 0.25*(solnData(VISC_VAR,i_ind-1,j_ind,k_ind-1) + &
                      solnData(VISC_VAR,i_ind-1,j_ind,k_ind)   + &
                      solnData(VISC_VAR,i_ind,j_ind,k_ind)     + &
                      solnData(VISC_VAR,i_ind,j_ind,k_ind-1))

           muwy = mu*zv(gridind)

           end select
           
        end if

        ! Release Pointer
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     enddo

  enddo

  ! Assign pressure and viscous forces to particleData:
  particleData(PRES_PART_PROP) = zL
  particleData(PEX0_PART_PROP) = particleData(PEXT_PART_PROP)
  particleData(PEXT_PART_PROP) = zpres

  ! Fvisc = -mu (N x W)
  particleData(FXVI_PART_PROP) = (nzp*muwy-nyp*muwz)
  particleData(FYVI_PART_PROP) = (nxp*muwz-nzp*muwx)
  particleData(FZVI_PART_PROP) = (nyp*muwx-nxp*muwy)

  ! Vorticity components:
  particleData(VORZ_PART_PROP) = zv(1)
  particleData(VORX_PART_PROP) = zv(2)
  particleData(VORY_PART_PROP) = zv(3)

  return

end subroutine ib_distributedForces


