!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/old_version2D/ib_forcing
!!
!!
!! NAME
!!
!! ib_forcing 
!!
!!
!! SYNOPSIS
!!
!!  ib_forcing_par(blockID, particleData,dt)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***
subroutine ib_forcing(blockID, particleData)


  ! This routine computes the immersed boundary forcing for a 
  ! given time step on a given particle.  A diffuse interface method is employed.  
  
  use Grid_data, ONLY: gr_meshMe
  
  use ImBound_data, only : ib_stencil,ib_alphax,ib_alphay
  
  use ib_interface, ONLY : ib_stencils,ib_interpLpoints,ib_extrapEpoints
  
  use Grid_interface, ONLY : Grid_getDeltas,          &
                             Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize, &
                             Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr

  use Driver_interface, ONLY : Driver_abortFlash

  use ImBound_Data, only : ib_dt
  
  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"
  
  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockID
  real, INTENT(INOUT) :: particleData(NPART_PROPS)
  !! -----------------------------------------------------
  
  integer, parameter :: ng = NGUARD
  
  integer, parameter :: nxc = NXB + NGUARD + 1
  integer, parameter :: nyc = NYB + NGUARD + 1
#if NDIM == 3
  integer, parameter :: nzc = NZB + NGUARD + 1
#else
  integer, parameter :: nzc = 1
#endif

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  ! counters 
  integer :: i,j,k
  
  real :: del(MDIM),coord(MDIM),bsize(MDIM)
  real :: invdt
  real :: ib_FuL,ib_FvL
  
  integer :: ierr

  real :: ib_dsxu,ib_dsyu,ib_dsxv,ib_dsyv,ib_dsxu2, ib_dsyv2  

  integer, dimension(ib_stencil) :: ib_ielemu,ib_jelemu,ib_ielemv, &
           ib_jelemv,ib_kelemu,ib_kelemv,ib_ielemw,ib_jelemw,ib_kelemw

  real, dimension(ib_stencil) :: ib_phileu,ib_philev,ib_philew

  real :: xb,yb,zb,sb,part_Velx,part_Vely,part_Velz,ib_ul,ib_vl,ib_wl

  ! Uncomment to set position and velocities to 0
  !spcoord = 0.0
  !ubdd = 0.0
  !vbdd = 0.0
  !ubd0 = 0.0
  !vbd0 = 0.0

  invdt = 1./ib_dt

  ! Indexes of stencils
  ib_ielemu = 0
  ib_jelemu = 0
  ib_ielemv = 0
  ib_jelemv = 0

#if NDIM == 3
  ib_kelemu = 0
  ib_kelemv = 0
  ib_ielemw = 0
  ib_jelemw = 0
  ib_kelemw = 0
#endif

  ib_dsxu=0.; ib_dsyu=0.
  ib_dsxv=0.; ib_dsyv=0.
  ib_dsxu2=0.; ib_dsyv2=0.


  ! Get face data (velocities):
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
  call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

  ! Get dx,dy
  call Grid_getDeltas(blockID,del)
  
  ! Get block center coords:
  call Grid_getBlkCenterCoords(blockID,coord)
  
  ! Get block size:
  call Grid_getBlkPhysicalSize(blockID,bsize)
!  
! Extract all the data from particleData

  xb = particleData(POSX_PART_PROP);
  yb = particleData(POSY_PART_PROP);
  sb = particleData(AREA_PART_PROP);
  part_Velx = particleData(VELX_PART_PROP);
  part_Vely = particleData(VELY_PART_PROP);
#if NDIM == 3
  zb = particleData(POSZ_PART_PROP);
  part_Velz = particleData(VELZ_PART_PROP);
#endif


#if NDIM == 2
  ! Obtain Stencils for interpolation to Lagrangian points:

  ! U velocities                                     
  if ((yb-coord(JAXIS)) .le. bsize(JAXIS)/2.+del(JAXIS)) then
  call ib_stencils(ng,nxc+1,nyc,xb,yb,ib_dsxu,ib_dsyu,    & 
                   ib_ielemu,ib_jelemu,IAXIS,             &
                   del,coord,bsize)

  ! Get the hL for that particle
  particleData(HL_PART_PROP)= (ib_dsxu+ib_dsyu)/(ib_alphax+ib_alphay)

  ! Interpolation of the values of velocity to Lagrangian points:
  call ib_interpLpoints(xb,yb,ib_dsxu,ib_dsyu,                       &
                        ib_ielemu,ib_jelemu,ib_phileu,ib_UL,          &
                        facexData(VELC_FACE_VAR,:,:,:),IAXIS,         &
                        GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,ng,nxc+1,nyc,del,coord,bsize)


  ! The particle forcing 
  ib_FuL  = invdt*(part_Velx - ib_UL ) 
  particleData(FUL_PART_PROP)  = particleData(FUL_PART_PROP) + ib_FuL

  ! Extrapolation of Forces to the Eulerian Points (Check what 
  ! happens to the positions accelerations and velocities) and 
  ! sum to ustar and wstar:
  call ib_extrapEpoints(xb,yb,sb,                                      &
                        ib_ielemu,ib_jelemu,ib_phileu,ib_FuL,          &
                        facexData(FORC_FACE_VAR,:,:,:),                &
                        GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,del,coord,bsize) 


  endif

  ! V velocities
  if ((xb-coord(IAXIS)) .le. bsize(IAXIS)/2.+del(IAXIS)) then
  call ib_stencils(ng,nxc,nyc+1,xb,yb,ib_dsxv,ib_dsyv,    &
                   ib_ielemv,ib_jelemv,JAXIS,             &
                   del,coord,bsize)

  ! Get the hL for that particle
  particleData(HL_PART_PROP)= (ib_dsxv+ib_dsyv)/(ib_alphax+ib_alphay)

  call ib_interpLpoints(xb,yb,ib_dsxv,ib_dsyv,                           &
                        ib_ielemv,ib_jelemv,ib_philev,ib_VL,             &
                        faceyData(VELC_FACE_VAR,:,:,:),JAXIS,            &
                        GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,ng,nxc,nyc+1,del,coord,bsize)                        

  ! The particle forcing
  ib_FvL = invdt*(part_Vely - ib_VL ) 
  particleData(FVL_PART_PROP)  = particleData(FVL_PART_PROP) + ib_FvL 
 
  ! Extrapolation of Forces to the Eulerian Points (Check what 
  ! happens to the positions accelerations and velocities) and 
  ! sum to ustar and wstar:
  call ib_extrapEpoints(xb,yb,sb,       &
                        ib_ielemv,ib_jelemv,ib_philev,ib_FvL,faceyData(FORC_FACE_VAR,:,:,:), &
                        GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,del,coord,bsize) 
                             
  endif

                                     
#elif NDIM == 3
  if (gr_meshMe .eq. MASTER_PE) write(*,*) 'Stencils not 3D yet!'
  call Driver_abortFlash("ib_forcing: Not 3D ready yet!")
#endif                                               
  
              
  ! Release face data (velocities):
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

end subroutine ib_forcing

