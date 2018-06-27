!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/rigid/sm_ioWriteStates_rigid
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!  Write IO variables and files for SolidMechanics rigid bodies.
!!
!!***

subroutine sm_ioWriteStates_rigid(ibd,nstep,time)

  use SolidMechanics_data, only: sm_bodyInfo, sm_structure

  implicit none
#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

  integer, intent(in) :: ibd,nstep
  real, intent(in)    :: time

  ! Internal variables:
  type(sm_structure),  pointer :: body

  real, dimension(NDIM) :: point, force_pres, force_visc, moment_pres, moment_visc

  integer :: imaster,i, idx, i_dim
  integer :: ix,ex,ia,ea,iw,ew

  real :: x_bod,y_bod,z_bod, xd_bod,yd_bod,zd_bod,xdd_bod,ydd_bod,zdd_bod

  character(len=6) :: str_ibd

  ! Get the body
  body => sm_BodyInfo(ibd)

  ! Extract info from Master:
  imaster = body % borigin_node 

  ix = body%ix; ex = body%ex;
  ia = body%ia; ea = body%ea;
  iw = body%iw; ew = body%ew;

  ! Dump forces into pres and visc variables:
  force_pres  = 0.; force_visc  = 0.;
  moment_pres = 0.; moment_visc = 0.;  
  
  ! Write to force_ and moment_ for writing to file out:
  ! Forces:
  i_dim = 0
  do i = body%ix,body%ex
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     force_pres(i_dim) = body%Hs_pres(idx); force_visc(i_dim) = body%Hs_visc(idx)
  end do

  ! Moments: 
  i_dim = 0
  do i = body%iw,body%ew
     i_dim = i_dim+1
     idx = body%ID(i,imaster)
     moment_pres(i_dim) = body%Hs_pres(idx); moment_visc(i_dim) = body%Hs_visc(idx)
  end do

  ! Translational variables:
  x_bod   = body%x(imaster)+body%qn(body%ID(ix,imaster))
  y_bod   = body%y(imaster)+body%qn(body%ID(ix+1,imaster))

  xd_bod  = body%qdn(body%ID(ix,imaster))
  yd_bod  = body%qdn(body%ID(ix+1,imaster))

  xdd_bod = body%qddn(body%ID(ix,imaster))
  ydd_bod = body%qddn(body%ID(ix+1,imaster))

#if NDIM == MDIM
  z_bod   = body%z(imaster)+body%qn(body%ID(ix+2,imaster))
  zd_bod  = body%qdn(body%ID(ix+2,imaster))
  zdd_bod = body%qddn(body%ID(ix+2,imaster))
#endif

  ! Write Out:
  write(str_ibd,"(I6.6)") ibd

  ! Write Forces, Moments:
#ifdef WRITEFORCES
  ! First call should done in SoliMechanics_init.F90
  open(unit=113,file='./IOData/force.'//str_ibd//'.res',form='formatted',&
       status='old',position='append')
  write(113,'(8g18.10)') nstep,time,force_pres(1:NDIM),force_visc(1:NDIM)
  close(113)
 
  open(unit=113,file='./IOData/momt.'//str_ibd//'.res',form='formatted',&
       status='old',position='append')
  write(113,'(8g18.10)') nstep,time,moment_pres(1:NDIM),moment_visc(1:NDIM)
  close(113)
#endif 

  ! Translational vars
#ifdef WRITESTATES
#if NDIM == MDIM
  open(unit=113,file='./IOData/posvelacc_x.'//str_ibd//'.res',form='formatted', &
       status='old',position='append')
  write(113,'(11g18.10)') nstep,time,  x_bod,  y_bod,  z_bod, &
                                      xd_bod, yd_bod, zd_bod, &
                                     xdd_bod,ydd_bod,zdd_bod
  close(113)

#else
  open(unit=113,file='./IOData/posvelacc_x.'//str_ibd//'.res',form='formatted', &
       status='old',position='append')
  write(113,'(11g18.10)') nstep,time,  x_bod,  y_bod, &
                                      xd_bod, yd_bod, &
                                     xdd_bod,ydd_bod
  close(113)

#endif
#endif

  ! Rotation variables
#ifdef WRITESTATES
  open(unit=113,file='./IOData/posvelacc_ang.'//str_ibd//'.res',form='formatted', &
       status='old',position='append')

  select case( body % trmatrix )
  case( RB_IDENTITY )
  case( RB_EULER321 )
  write(113,'(11g18.10)') nstep,time,  body%qn(body%ID(ia:ea,imaster)), &
                                      body%qdn(body%ID(ia:ea,imaster)), &
                                     body%qddn(body%ID(ia:ea,imaster))
  case( RB_QUATERNN )
     call Driver_abortFlash("sm_ioWriteStates_rigid: Quaternion write not implemented")
  case( RB_TWODIM )
  write(113,'(11g18.10)') nstep,time,  body%qn(body%ID(ia:ea,imaster)), &
                                      body%qdn(body%ID(ia:ea,imaster)), &
                                     body%qddn(body%ID(ia:ea,imaster))
  case default
     call Driver_abortFlash("sm_ioWriteStates_rigid: trmatrix Type not known")
  end select
  close(113)
#endif


  ! Angular velocities
#ifdef WRITESTATES
  open(unit=113,file='./IOData/posvelacc_w.'//str_ibd//'.res',form='formatted', &
       status='old',position='append')

  select case( body % trmatrix )
  case( RB_EULER321 , RB_IDENTITY )
  write(113,'(11g18.10)') nstep,time,  body%qn(body%ID(iw:ew,imaster)), &
                                      body%qdn(body%ID(iw:ew,imaster)), &
                                     body%qddn(body%ID(iw:ew,imaster))
  case( RB_QUATERNN )
     call Driver_abortFlash("sm_ioWriteStates_rigid: Quaternion write not implemented")
  case( RB_TWODIM )
  write(113,'(11g18.10)') nstep,time,  body%qn(body%ID(iw:ew,imaster)), &
                                      body%qdn(body%ID(iw:ew,imaster)), &
                                     body%qddn(body%ID(iw:ew,imaster))
  case default
     call Driver_abortFlash("sm_ioWriteStates_rigid: trmatrix Type not known")
  end select
  close(113)
#endif

  return

end subroutine sm_ioWriteStates_rigid
