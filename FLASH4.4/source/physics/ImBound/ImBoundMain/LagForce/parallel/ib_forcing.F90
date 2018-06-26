!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ib_forcing
!!
!!
!! NAME
!!
!! ib_forcing 
!!
!!
!! SYNOPSIS
!!
!!  ib_forcing(blockID, particleData)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***
subroutine ib_forcing(ibd,p,blockID,particleData)

  ! This routine computes the immersed boundary forcing for a 
  ! given time step on a given particle.  A diffuse interface method is employed.
  use Grid_data, ONLY: gr_meshMe
  
  use ImBound_data, only : ib_stencil,ib_alphax,ib_alphay,ib_dt,ib_BlockMarker_flag
  
  use ib_interface, ONLY : ib_stencils,      &
                           ib_interpLpoints, &
                           ib_extrapEpoints
  
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
   
  integer,parameter,dimension(MDIM):: FORCE_IND=(/FUL_PART_PROP,FVL_PART_PROP,FWL_PART_PROP/)
  integer,parameter,dimension(MDIM):: VELITP_IND=(/UITP_PART_PROP,VITP_PART_PROP,WITP_PART_PROP/) ! Shizhao Wang
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

  invdt = 1./ib_dt

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
  part_Vel(IAXIS) = particleData(VELX_PART_PROP);
  part_Vel(JAXIS) = particleData(VELY_PART_PROP);

  part_Nml(IAXIS) = particleData(NMLX_PART_PROP);
  part_Nml(JAXIS) = particleData(NMLY_PART_PROP);

#if NDIM == 3
  part_Pos(KAXIS) = particleData(POSZ_PART_PROP);
  part_Vel(KAXIS) = particleData(VELZ_PART_PROP);
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
        ib_FuL(i)  = invdt*(part_Vel(i) - ib_UL(i)) 
        particleData(FORCE_IND(i))  = particleData(FORCE_IND(i)) + ib_FuL(i)
        particleData(VELITP_IND(i))  = ib_uL(i)

        particleData(HL_PART_PROP) = hl

        !call Timers_start("ib_extrapEpoints")
        ! Extrapolation of Forces to the Eulerian Points (Check what 
        ! happens to the positions accelerations and velocities) and 
        ! sum to ustar, vstar and wstar:
        call ib_extrapEpoints(part_Pos,sb,hl,del,ib_ielem(:,:,i),ib_phile(:,:,i),   &        
                              ib_FuL(i),blockID,FACE_IND(i))
        !call Timers_stop("ib_extrapEpoints")

     endif

  end do

  elseif(gr_sbFirstcall .eq. CONSTANT_ONE) then ! Fixed bod and First call

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
        ib_FuL(i)  = invdt*(part_Vel(i) - ib_UL(i))
        particleData(FORCE_IND(i))  = particleData(FORCE_IND(i)) + ib_FuL(i)

        particleData(HL_PART_PROP) = hl

        !call Timers_start("ib_extrapEpoints")
        ! Extrapolation of Forces to the Eulerian Points (Check what 
        ! happens to the positions accelerations and velocities) and 
        ! sum to ustar, vstar and wstar:
        call ib_extrapEpoints(part_Pos,sb,hl,del,ib_ielem(:,:,i),ib_phile(:,:,i),   &
                              ib_FuL(i),blockID,FACE_IND(i))
        !call Timers_stop("ib_extrapEpoints")

     endif

  end do

  ! Assign ielem and phile to gr_sbBodyInfo
  gr_sbBodyInfo(ibd)%ielem(:,:,:,p) = ib_ielem(:,:,:)
  gr_sbBodyInfo(ibd)%phile(:,:,:,p) = ib_phile(:,:,:)


  else ! Fixed bod and ~ first call

  ib_ielem(:,:,:) = gr_sbBodyInfo(ibd)%ielem(:,:,:,p)
  ib_phile(:,:,:) = gr_sbBodyInfo(ibd)%phile(:,:,:,p)

  do idir = 1,NDIM

    if (force_fl(idir)) then

    ! Get Pointer to faceData in blockID, direction X,Y,Z:
    call Grid_getBlkPtr(blockID,faceData,FACE_IND(idir))

    ! Interpolate
    ! Value of the function in xp,yp,zp:
    ib_uL(idir) = 0.;
    do i = 1 , ib_stencil
       ib_uL(idir) = ib_uL(idir) + ib_phile(i,CONSTANT_ONE,idir) * &
                     faceData(VELC_FACE_VAR,ib_ielem(i,IAXIS,idir),ib_ielem(i,JAXIS,idir),ib_ielem(i,KAXIS,idir));
    enddo
  
    ! Particle Forcing field:
    ib_FuL(idir)  = invdt*(part_Vel(idir) - ib_UL(idir))
    particleData(FORCE_IND(idir))  = particleData(FORCE_IND(idir)) + ib_FuL(idir)

    ! Extrapolation of Forces to the Eulerian Points
    factor = sb*particleData(HL_PART_PROP)/PRODUCT(del(1:NDIM))

    ! Do Extrapolation:
    do i = 1,ib_stencil

     faceData(FORC_FACE_VAR,ib_ielem(i,IAXIS,idir),ib_ielem(i,JAXIS,idir),ib_ielem(i,KAXIS,idir)) = &
     faceData(FORC_FACE_VAR,ib_ielem(i,IAXIS,idir),ib_ielem(i,JAXIS,idir),ib_ielem(i,KAXIS,idir)) + &
     factor*ib_phile(i,CONSTANT_ONE,idir)*ib_FuL(idir);

    enddo

    ! Release face data:
    call Grid_releaseBlkPtr(blockID,faceData,FACE_IND(idir))
    
    endif

  enddo

  !write(filename,"(A,I5.5)") "./forc", gr_meshMe
  !if (firstcall) then
  !  open(33,file=trim(filename),status='unknown')
  !  firstcall = .false.
  !else
  !  open(33,file=trim(filename), status='old', position='append')
  !endif
  !write(33,*) ibd,p,ib_uL(IAXIS:KAXIS),ib_FuL(IAXIS:KAXIS)
  !close(33)

  endif

  return

end subroutine ib_forcing




