!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/gr_sbDistributedForces
!!
!! NAME
!!  gr_sbDistributedForces
!!
!! SYNOPSIS
!!
!!  gr_sbDistributedForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbDistributedForces()
#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs, gr_maxParticlesPerProc
#else
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
#endif

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, &
                             Grid_getDeltas


  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbParticleCount, solid_body

  use gr_sbInterface, ONLY : gr_sbSendForces

  use ib_interface, ONLY : ib_distributedForces

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  type(solid_body), pointer :: bodyInfo
  integer, save, dimension(MAXBLOCKS) :: listOfBlocks
  integer, save :: count
  integer :: i,j,p,b,lb,blockID,localPart,gettingFrom

  real :: particleData(NPART_PROPS)

  real, dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1,MAXBLOCKS) :: vortz
  real, dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1,MAXBLOCKS) :: vortx
  real, dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1,MAXBLOCKS) :: vorty

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,solnData,facezData

  real :: del(MDIM),dx,dy,dz

  integer :: ix,ex,iy,ey,iz,ez


  ! Compute vorticity 
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! Get dx,dy,dz
     call Grid_getDeltas(blockID,del)
                 
     ! Calculate the vorticity of the block:
     dx = del(IAXIS) 
     dy = del(JAXIS) 
     dz = 1.
#if NDIM == MDIM
     dz = del(KAXIS)
#endif

     ! Initialize vorticity:
#if NDIM == MDIM
     vortx(:,:,:,blockID) = 0.
     vorty(:,:,:,blockID) = 0.
#endif
     vortz(:,:,:,blockID) = 0.

     ! Compute Vorticities: Works only for two or more layers of Guardcells
     ! Compute wz
     ix = (CONSTANT_TWO-1)*K1D  + 1
     ex = (GRID_IHI_GC-1) *K1D  + 1
     iy = (CONSTANT_TWO-1)*K2D  + 1
     ey = (GRID_JHI_GC-1) *K2D  + 1
     iz =  CONSTANT_ONE 
     ez = (GRID_KHI_GC-1) *K3D  + 1

     vortz(ix:ex,iy:ey,iz:ez,blockID) = &                             ! dv/dx-du/dy
         ( faceyData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,iz:ez) -       &
           faceyData(VELC_FACE_VAR,ix-1:ex-1,iy:ey    ,iz:ez) )/dx -  &
         ( facexData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,iz:ez) -       &
           facexData(VELC_FACE_VAR,ix:ex    ,iy-1:ey-1,iz:ez) )/dy
     vortz(CONSTANT_ONE ,iy:ey    ,iz:ez,blockID) = vortz(CONSTANT_TWO,iy:ey    ,iz:ez,blockID)
     vortz(GRID_IHI_GC+1,iy:ey    ,iz:ez,blockID) = vortz(GRID_IHI_GC ,iy:ey    ,iz:ez,blockID)
     vortz(ix:ex  ,CONSTANT_ONE   ,iz:ez,blockID) = vortz(ix:ex ,CONSTANT_TWO   ,iz:ez,blockID)
     vortz(ix:ex  ,GRID_JHI_GC+1  ,iz:ez,blockID) = vortz(ix:ex ,GRID_JHI_GC    ,iz:ez,blockID)


#if NDIM == MDIM
     ! Compute wx, wy:
     ix = CONSTANT_ONE
     iz = (CONSTANT_TWO-1)*K3D  + 1
     vortx(ix:ex,iy:ey,iz:ez,blockID) = &                                 ! dw/dy-dv/dz
         ( facezData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,    iz:ez)      -  &
           facezData(VELC_FACE_VAR,ix:ex    ,iy-1:ey-1,    iz:ez) )/dy -  &
         ( faceyData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,    iz:ez) -       &
           faceyData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,iz-1:ez-1) )/dz
     vortx(ix:ex  ,CONSTANT_ONE   ,iz:ez,blockID) = vortx(ix:ex ,CONSTANT_TWO   ,iz:ez,blockID)
     vortx(ix:ex  ,GRID_JHI_GC+1  ,iz:ez,blockID) = vortx(ix:ex ,GRID_JHI_GC    ,iz:ez,blockID)
     vortx(ix:ex  ,iy:ey    ,CONSTANT_ONE   ,blockID) = vortx(ix:ex ,iy:ey ,CONSTANT_TWO   ,blockID)
     vortx(ix:ex  ,iy:ey    ,GRID_KHI_GC+1  ,blockID) = vortx(ix:ex ,iy:ey ,GRID_KHI_GC    ,blockID)


     iy = CONSTANT_ONE
     ix = (CONSTANT_TWO-1)*K1D  + 1
     vorty(ix:ex,iy:ey,iz:ez,blockID) = &                                 ! du/dz-dw/dx
         ( facexData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,    iz:ez)      -  &
           facexData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,iz-1:ez-1) )/dz -  &
         ( facezData(VELC_FACE_VAR,ix:ex    ,iy:ey    ,    iz:ez) -       &
           facezData(VELC_FACE_VAR,ix-1:ex-1,iy:ey    ,    iz:ez) )/dx
     vorty(CONSTANT_ONE ,iy:ey    ,iz:ez,blockID) = vorty(CONSTANT_TWO,iy:ey    ,iz:ez,blockID)
     vorty(GRID_IHI_GC+1,iy:ey    ,iz:ez,blockID) = vorty(GRID_IHI_GC ,iy:ey    ,iz:ez,blockID)
     vorty(ix:ex  ,iy:ey    ,CONSTANT_ONE   ,blockID) = vorty(ix:ex ,iy:ey ,CONSTANT_TWO   ,blockID)
     vorty(ix:ex  ,iy:ey    ,GRID_KHI_GC+1  ,blockID) = vorty(ix:ex ,iy:ey ,GRID_KHI_GC    ,blockID)
#endif

     ! Get face data (velocities):
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo



  ! Loop over all bodies to get pressure and viscous forces in markers:
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     ! Master Proc:
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        ! Find number of particles in the Master Processor:
        localPart = 0    
        do j=1,size(bodyInfo%particlesPerProc,DIM=2)
           if(bodyInfo%particlesPerProc(1,j) .eq. bodyInfo%bodyMaster) then
              localPart = bodyInfo%particlesPerProc(2,j)
           endif
        enddo

        ! Initialize Particle distributed Forces:
        bodyInfo % particles(PRES_PART_PROP,:) = 0.
        bodyInfo % particles(FXVI_PART_PROP,:) = 0.
        bodyInfo % particles(FYVI_PART_PROP,:) = 0.
#if NDIM == 3
        bodyInfo % particles(FZVI_PART_PROP,:) = 0.
#endif

        do i = 1, localPart
           if (int(bodyInfo%particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then

              blockID = int(bodyInfo%particles(BLK_PART_PROP,i))
              particleData = bodyInfo%particles(1:NPART_PROPS,i)

              call ib_distributedForces(blockID, particleData, vortx(:,:,:,blockID), vorty(:,:,:,blockID), vortz(:,:,:,blockID))

              bodyInfo % particles(1:NPART_PROPS,i) = particleData   

           end if
        end do

     ! All other Procs:
     else

        gettingFrom = gr_sbParticleCount(b)

        if (gettingFrom .gt. 0) then

           ! Initialize Particle distributed Forces:
           bodyInfo % particles(PRES_PART_PROP,:) = 0.
           bodyInfo % particles(FXVI_PART_PROP,:) = 0.
           bodyInfo % particles(FYVI_PART_PROP,:) = 0.
#if NDIM == 3
           bodyInfo % particles(FZVI_PART_PROP,:) = 0.
#endif

           do p = 1,gettingFrom

              blockID = int(bodyInfo%particles(BLK_PART_PROP,p))
              particleData = bodyInfo%particles(1:NPART_PROPS,p)

              call ib_distributedForces(blockID, particleData, vortx(:,:,:,blockID), vorty(:,:,:,blockID), vortz(:,:,:,blockID))

              bodyInfo % particles(1:NPART_PROPS,p) = particleData 

           enddo
        endif
     endif
     nullify(bodyInfo)
  enddo



  ! Now gather all the particle info back to the Masters so they
  ! can do the integration of distributed forces.
  call gr_sbSendForces()



  
!!$  ! Loop over all bodies
!!$  do b = 1, gr_sbNumBodies
!!$     bodyInfo => gr_sbBodyInfo(b)
!!$     ! Master Proc:
!!$     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!!$        open(33,file="./force2.dat",status="replace")
!!$        do i = 1, bodyInfo%totalPart, 2    
!!$           write(33,'(2I4,3g18.10,I4)') i,bodyInfo % bodyMaster, &
!!$                      bodyInfo%particles(PRES_PART_PROP,i),bodyInfo%particles(FXVI_PART_PROP,i),bodyInfo%particles(FYVI_PART_PROP,i),int(bodyInfo%particles(GLOB_PART_PROP,i))
!!$        enddo
!!$        close(33)
!!$     endif
!!$     nullify(bodyInfo)
!!$  enddo

  return

end subroutine gr_sbDistributedForces
