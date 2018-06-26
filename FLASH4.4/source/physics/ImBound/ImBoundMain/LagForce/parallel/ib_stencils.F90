!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_stencils
!!
!!
!! NAME
!!
!! ib_stencils 
!!
!!
!! SYNOPSIS
!!
!!  subroutine ib_stencils(xp,np,gridfl,del,coord,bsize,   & 
!!                         ielem,hl,forceflag)
!!                          
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ib_stencils(xp,np,gridfl,del,coord,bsize,   & 
                       ielem,hl,forcflag)


  use ImBound_data , ONLY : ib_interp, ib_stencil, ib_alphax, ib_alphay, ib_alphaz

  implicit none
#include "Flash.h"
#include "constants.h"
#include "ImBound.h"
  
  integer, INTENT(IN) :: gridfl(MDIM),forcflag
  real, INTENT(IN) :: xp(MDIM),np(MDIM)
  integer, INTENT(INOUT) :: ielem(ib_stencil,MDIM)
  real, INTENT(OUT) :: hl
  real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
   
  ! Local Variables:
  integer :: i,idim

  integer, parameter :: ng = NGUARD
  integer, parameter :: nxc = NXB + NGUARD + 1
  integer, parameter :: nyc = NYB + NGUARD + 1
#if NDIM == 3
  integer, parameter :: nzc = NZB + NGUARD + 1
#else
  integer, parameter :: nzc = 1
#endif
  integer, parameter, dimension(MDIM) :: nc = (/nxc,nyc,nzc/)
  integer, parameter, dimension(MDIM) :: nf = (/GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC/)
  
  real    :: auxx,dist(MDIM)
  integer :: cpoint(MDIM),icpoint,jcpoint,kcpoint
  integer :: bnd(LOW:HIGH,MDIM)
  real,save :: dsx,dsy,dsz,hli
  real    :: delaux(MDIM),xi
  real    :: eps

#if NDIM == 2 
  integer, parameter :: indIaxis(1:9) = (/ 0, 0, 0, 1, -1, 1, -1, -1, 1 /)
  integer, parameter :: indJaxis(1:9) = (/ 0, 1,-1, 0,  0, 1,  1, -1,-1 /)
#elif NDIM ==3 
  integer, parameter :: indIaxis(1:15) = (/ 0, 1,-1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1 /)
  integer, parameter :: indJaxis(1:15) = (/ 0, 0, 0, 1,-1, 0, 0, 1,  1, -1,-1, 1,  1, -1,-1 /)
  integer, parameter :: indKaxis(1:15) = (/ 0, 0, 0, 0, 0, 1,-1, 1,  1,  1, 1,-1, -1, -1,-1 /)
#endif

  integer :: ipt
  real :: loc_start, loc_factor
  logical, save :: firstcall=.true.


  cpoint(:) = 0

#ifdef IB_OPTIMIZE


  ! Find closest point:
  ! Closest point in x,y,z
  do idim = 1,NDIM

     loc_start    = coord(idim) - 0.5*bsize(idim)
     loc_factor   = (xp(idim) - loc_start)/del(idim)

     ! Use max CONSTANT_TWO to test virtual particles low side:
     cpoint(idim) = max(CONSTANT_TWO,(floor(loc_factor) + NGUARD + 1))

     if (gridfl(idim) .eq. FACES) then

        xi= loc_start + real(cpoint(idim) - NGUARD - 1)*del(idim)

        ! Virtual particle test:
        if (cpoint(idim) .eq. nf(idim)) cycle  

        if ( abs(xp(idim)-xi) .gt. (0.5*del(idim)) ) cpoint(idim) = cpoint(idim)+1  
        
     else
        ! Virtual particle test:
        !if (cpoint(idim) .gt. nc(idim)) cpoint(idim) = nc(idim)
        if (cpoint(idim) .gt. (nf(idim)-1)) cpoint(idim) = nf(idim)-1 ! This is s.t. when using CF for
                                                                      ! shear stress computation on IB
                                                                      ! The range is extended up to the
                                                                      ! previous from last guardcell point.
 
     endif
  enddo

#else

  ! Find closest point: Old Beastly Way
  !eps = 1.E-2*MINVAL(del(1:NDIM))
  bnd(:,:)  =   0
  select case (forcflag)
  case(FORCE_FLOW)

     do i=1,NDIM

        bnd(LOW,i)  = NGUARD
        bnd(HIGH,i) = nc(i)
        delaux(i)   = 0.5*del(i)
        if (gridfl(i) .eq. FACES) then
           bnd(HIGH,i) = bnd(HIGH,i) + 1
           delaux(i) = 0.
        endif

     enddo

  case(COMPUTE_FORCES)

     do i=1,NDIM
        
        bnd(LOW,i)  = CONSTANT_TWO
        bnd(HIGH,i) = nf(i)
        delaux(i) = 0.
        if (gridfl(i) .eq. CENTER) then
           bnd(HIGH,i) = bnd(HIGH,i)- 1 ! Pressure
           delaux(i)   = 0.5*del(i)
        endif
     enddo

  case default

     call Driver_abortFlash("ib_sencils : forcflag does not correspond to any available option.") 
  
  end select

  dist(:)   = 10.e10
  do idim = 1,NDIM
     do i = bnd(LOW,idim),bnd(HIGH,idim)
        xi = coord(idim) - 0.5*bsize(idim) + &
             real(i - ng - 1)*del(idim) + delaux(idim)
        auxx = abs(xp(idim) - xi);
        if (auxx .lt. dist(idim)) then !+eps Took out this epsilon addition as biases stencil 
                                       !     location for marker right in the middle between
                                       !     two velocity points. 
           cpoint(idim) = i;
           dist(idim) = auxx;
        endif
     enddo
  enddo 
 
#endif

  
  ! dsx, dsy, dsz:
#ifdef FLASH_GRID_PARAMESH
  dsx = ib_alphax*del(IAXIS)
  dsy = ib_alphay*del(JAXIS)
#if NDIM == 2
  hl  = (dsx+dsy)/(ib_alphax+ib_alphay)
#elif NDIM == 3
  dsz = ib_alphaz*del(KAXIS)
  hl  = (dsx+dsy+dsz)/(ib_alphax+ib_alphay+ib_alphaz)
#endif
#else
  if (firstcall) then
     dsx = ib_alphax*del(IAXIS)
     dsy = ib_alphay*del(JAXIS)
#if NDIM == 2
     !hli  = (dsx+dsy)/(ib_alphax+ib_alphay)
#elif NDIM == 3
     dsz = ib_alphaz*del(KAXIS)
     !hli  = (dsx+dsy+dsz)/(ib_alphax+ib_alphay+ib_alphaz)
#endif
     firstcall = .false.
  endif
  !hl = hli
  ! Case of distorted cells in UG mode: 
#if NDIM == 2
  hl = sqrt((del(IAXIS)*np(IAXIS))**2. + (del(JAXIS)*np(JAXIS))**2.)
#elif NDIM == 3
  hl = sqrt((del(IAXIS)*np(IAXIS))**2. + (del(JAXIS)*np(JAXIS))**2. + (del(KAXIS)*np(KAXIS))**2.)
#endif

#endif


  ! Build stencil structure:
  ielem(1:ib_stencil,1:MDIM) = 1

#if NDIM == 2

  do ipt=1,ib_stencil

     ielem(ipt,IAXIS) = cpoint(IAXIS) + indIaxis(ipt)
     ielem(ipt,JAXIS) = cpoint(JAXIS) + indJaxis(ipt)

  enddo

!!$  select case (ib_interp)
!!$  case(CONSTANT_ONE) ! Linear basis polynomial.
!!$     ! Five Point Stencil
!!$     ! Center North South East West:
!!$     ielem(1:ib_stencil,IAXIS) = (/ icpoint,icpoint,icpoint,     &
!!$                                    icpoint+1,icpoint-1 /)
!!$     ielem(1:ib_stencil,JAXIS) = (/ jcpoint,jcpoint+1,jcpoint-1, &
!!$                                    jcpoint,jcpoint /)
!!$ 
!!$  case(CONSTANT_TWO) ! Quadratic basis polynomial.
!!$     ! Nine Point Stencil
!!$     ! Center North South East West NorthEast NorthWest SouthWest SouthEast
!!$     if (ib_stencil .eq. 9) then
!!$     ielem(1:ib_stencil,IAXIS) = (/ icpoint,icpoint,icpoint,        &
!!$                                    icpoint+1,icpoint-1, icpoint+1, &
!!$                                    icpoint-1, icpoint-1, icpoint+1 /)
!!$     ielem(1:ib_stencil,JAXIS) = (/ jcpoint,jcpoint+1,jcpoint-1, &
!!$                                    jcpoint,jcpoint, jcpoint+1,  &
!!$                                    jcpoint+1, jcpoint-1, jcpoint-1 /)     
!!$     endif
!!$
!!$  end select

#elif NDIM == 3


  do ipt=1,ib_stencil

     ielem(ipt,IAXIS) = cpoint(IAXIS) + indIaxis(ipt)
     ielem(ipt,JAXIS) = cpoint(JAXIS) + indJaxis(ipt)
     ielem(ipt,KAXIS) = cpoint(KAXIS) + indKaxis(ipt)

  enddo


!!$  select case(ib_interp)
!!$  case(CONSTANT_ONE) 
!!$     ! 7 point stencil
!!$     ! Center East West North South Back Front:    
!!$     ielem(1:ib_stencil,IAXIS)=(/ icpoint,icpoint+1,icpoint-1,   &
!!$                                  icpoint,icpoint,icpoint,icpoint /)
!!$
!!$     ielem(1:ib_stencil,JAXIS) = (/ jcpoint,jcpoint,jcpoint,     &
!!$                                    jcpoint+1,jcpoint-1,jcpoint,jcpoint /)
!!$
!!$     ielem(1:ib_stencil,KAXIS) = (/ kcpoint,kcpoint,kcpoint,     &
!!$                                    kcpoint,kcpoint,kcpoint+1,kcpoint-1 /)  
!!$  case(CONSTANT_TWO)
!!$     ! 15 point stencil:
!!$     ielem(1:ib_stencil,IAXIS) = (/ icpoint,icpoint+1,icpoint-1,               &
!!$                                    icpoint,icpoint,icpoint,icpoint,           &
!!$                                    icpoint+1,icpoint-1,icpoint-1,icpoint+1,   &
!!$                                    icpoint+1,icpoint-1,icpoint-1,icpoint+1 /)
!!$     ielem(1:ib_stencil,JAXIS) = (/ jcpoint,jcpoint,jcpoint,                   &
!!$                                    jcpoint+1,jcpoint-1,jcpoint,jcpoint,       &
!!$                                    jcpoint+1,jcpoint+1,jcpoint-1,jcpoint-1,   &
!!$                                    jcpoint+1,jcpoint+1,jcpoint-1,jcpoint-1 /)
!!$     ielem(1:ib_stencil,KAXIS) = (/ kcpoint,kcpoint,kcpoint,                   &
!!$                                    kcpoint,kcpoint,kcpoint+1,kcpoint-1,       &
!!$                                    kcpoint+1,kcpoint+1,kcpoint+1,kcpoint+1,   &
!!$                                    kcpoint-1,kcpoint-1,kcpoint-1,kcpoint-1 /)
!!$  
!!$  end select



#endif 
  

  return

end subroutine ib_stencils







