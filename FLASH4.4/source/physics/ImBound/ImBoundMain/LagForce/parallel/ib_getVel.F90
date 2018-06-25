!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_getVel
!!
!!
!! NAME
!!
!! ib_getVel
!!
!!
!! SYNOPSIS
!!
!!  ib_getVel(blockID, particleData)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***
subroutine ib_getVel(ibd,p,blockID,particleData)

  ! This routine interpolates the velocity on Eulerian points to Lagrangian point  
  ! given time step on a given particle.  A diffuse interface method is employed.
  use Grid_data, ONLY: gr_meshMe
  
  use ImBound_data, only : ib_stencil,ib_alphax,ib_alphay,ib_dt,ib_BlockMarker_flag
  
  use ib_interface, ONLY : ib_stencils,      &
                           ib_interpLpoints
  
  use Grid_interface, ONLY : Grid_getDeltas,          &
                             Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize, &
                             Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr  

  use Driver_interface, ONLY : Driver_abortFlash

  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  use gr_sbData, only : gr_sbBodyInfo, gr_sbFirstCall

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"
#include "ImBound.h"
  
  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockID,ibd,p
  real, INTENT(INOUT) :: particleData(NPART_PROPS)
  !! -----------------------------------------------------
   
!  integer,parameter,dimension(MDIM):: FORCE_IND=(/FUL_PART_PROP,FVL_PART_PROP,FWL_PART_PROP/)
   !integer,parameter,dimension(MDIM):: VELITP_IND=(/UITP_PART_PROP,VITP_PART_PROP,WITP_PART_PROP/) ! Shizhao Wang
   integer,parameter,dimension(MDIM):: VELITP_IND=(/UUPD_PART_PROP,VUPD_PART_PROP,WUPD_PART_PROP/) ! Shizhao Wang
   integer,parameter,dimension(MDIM):: FACE_IND =(/FACEX,FACEY,FACEZ/)

  ! counters 
  integer :: i,j,k,idir
  
  real :: del(MDIM),coord(MDIM),bsize(MDIM)
  real :: invdt
  
  integer :: ierr

  integer, dimension(ib_stencil,MDIM,MDIM) :: ib_ielem

  real, dimension(ib_stencil,NDIM+1,MDIM) :: ib_phile

  real :: part_Pos(MDIM),part_Vel(MDIM),sb,hl,ib_uL(MDIM),ib_FuL(MDIM)

  real :: part_Nml(MDIM)

  integer :: indx(MDIM),gridfl(MDIM)

  logical :: force_fl(MDIM),force_dir(MDIM)

  real, pointer, dimension(:,:,:,:) :: faceData

  real :: factor

  !character(len=20) :: filename
  !logical, save :: firstcall = .true.

  ib_BlockMarker_flag(blockID) = .TRUE.

  ! Get dx,dy,dz
  call Grid_getDeltas(blockID,del)
  
  ! Get block center coords:
  call Grid_getBlkCenterCoords(blockID,coord)
  
  ! Get block size:
  call Grid_getBlkPhysicalSize(blockID,bsize)

  ! Extract all the data from particleData
  part_Pos(:) = 0.
  part_Vel(:) = 0.
  part_Nml(:) = 0.
  sb = 0.

  part_Pos(IAXIS) = particleData(POSX_PART_PROP);
  part_Pos(JAXIS) = particleData(POSY_PART_PROP);
  sb = particleData(AREA_PART_PROP);
!  part_Vel(IAXIS) = particleData(VELX_PART_PROP);
!  part_Vel(JAXIS) = particleData(VELY_PART_PROP);

  part_Nml(IAXIS) = particleData(NMLX_PART_PROP);
  part_Nml(JAXIS) = particleData(NMLY_PART_PROP);

#if NDIM == 3
  part_Pos(KAXIS) = particleData(POSZ_PART_PROP);
!  part_Vel(KAXIS) = particleData(VELZ_PART_PROP);
  part_Nml(KAXIS) = particleData(NMLZ_PART_PROP);
#endif


  ! Check if particle is within the bounding box of each velocity grid
  force_fl(1:MDIM) = .true.
  force_dir(1:MDIM)= .true. 
  do i = 1,NDIM
     force_dir(i) =(part_Pos(i)-coord(i)) .le.  (bsize(i)/2. + del(i))  
  enddo
  force_fl(IAXIS) = force_dir(JAXIS) .and. force_dir(KAXIS)
  force_fl(JAXIS) = force_dir(IAXIS) .and. force_dir(KAXIS)
#if NDIM == 3
  force_fl(KAXIS) = force_dir(IAXIS) .and. force_dir(JAXIS)
#endif

  if (gr_sbBodyInfo(ibd)%sbIsFixed .eq. CONSTANT_ZERO) then
 
 
  do i = 1,NDIM
     ! Obtain Stencils for interpolation to Lagrangian points:
     ! U, V and W velocities
     gridfl(:) = CENTER
     indx(:)   = CONSTANT_ZERO      
                               
     if (force_fl(i)) then

        gridfl(i) = FACES
        indx(i) = CONSTANT_ONE

        !call Timers_start("ib_stencils")
        ! Define Interpolation Stencil For Particle:
        call ib_stencils(part_Pos,part_Nml,gridfl,del,coord,bsize,   & 
                         ib_ielem(:,:,i),hl,FORCE_FLOW)
        !call Timers_stop("ib_stencils")

        !call Timers_start("ib_interpLpoints")
        ! Interpolation of the values of velocity to Lagrangian points:
        call ib_interpLpoints(part_Pos,gridfl,                     &
              del,coord,bsize,ib_ielem(:,:,i),ib_phile(:,:,i),     &
              ib_uL(i),FORCE_FLOW,blockID,FACE_IND(i))
        !call Timers_stop("ib_interpLpoints")

        ! The particle forcing field:
        particleData(VELITP_IND(i))  = ib_uL(i)

     endif

  end do

  else ! Fixed bod and ~ first call
    write(*,*) 'Undefined paramters for ib_getVel'
    pause
  endif

  return

end subroutine ib_getVel




