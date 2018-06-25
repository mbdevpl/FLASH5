!!****if* source/physics/ImBound/ImBoundMain/LagForce/serial/ImBound_data
!!
!! NAME
!!
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data specific to the Immersed Boundary module.
!!
!!***


module ImBound_data

#include "Flash.h"
#include "constants.h"

      ! Number of bodies
      integer, save :: ib_nbd

      ! Maximum number of marker points in the system
      integer, parameter :: ib_nmaxa = 3000
      integer, save :: ib_nnoda

#if NDIM == 3
      ! Maximum number of aero elements that reach the node.
      integer, parameter :: ib_nmaxelnod = 15 + 1
#endif

      ! Aerodinamic surfaces elements
      integer, parameter :: ib_emaxa = 3000
      integer, parameter :: ib_nodelema = 4 ! Nodes per aerodynamic element + 1
      integer, save :: ib_nela

      ! Arrays with coordinates of marker points:
      real, save :: ib_xbus(ib_nmaxa),ib_ybus(ib_nmaxa),ib_xbcomp(ib_nmaxa),ib_ybcomp(ib_nmaxa), &
                    ib_xb(ib_nmaxa),ib_yb(ib_nmaxa)
      ! Arrays with marker points velocities and accelerations:
      real, save :: ib_ubd(ib_nmaxa),ib_vbd(ib_nmaxa),ib_ubdd(ib_nmaxa),ib_vbdd(ib_nmaxa),       &
                    ib_ubd0(ib_nmaxa),ib_vbd0(ib_nmaxa)

      ! Normals to the Lagrangian points:
      real, save :: ib_nxL(ib_nmaxa),ib_nyL(ib_nmaxa)
    
      ! Integers that define if a given grid block has markers inside
      integer, save :: ib_bflags(MDIM,MAXBLOCKS)  

      !
      integer, dimension(ib_nmaxa) :: ib_lpindexu,ib_lpindexv


#if NDIM == 3
      real, save :: ib_zbus(ib_nmaxa),ib_zbcomp(ib_nmaxa),ib_zb(ib_nmaxa)
      real, save :: ib_wbd(ib_nmaxa),ib_wbdd(ib_nmaxa),ib_wbd0(ib_nmaxa)
      real, save :: ib_nzL(ib_nmaxa)

      ! Aerodynamic elements normals:
      real, save :: ib_nexL(ib_emaxa),ib_neyL(ib_emaxa),ib_nezL(ib_emaxa)

      ! Array of aerodynamic elements conected to the aerodynamic nodes:
      integer, save :: ib_AELEMNODE(ib_nmaxelnod,ib_nmaxa)          ! PROBABLY DON'T NEED THIS

      ! Array of node associated areas:
      real, save :: ib_AAREANODE(ib_nmaxa),ib_AANGNODE(ib_nmaxa)

      ! Aerodynamic elements corner angles:
      real, save :: ib_AANGELEM(ib_nodelema-1,ib_emaxa)

      !
      integer, dimension(ib_nmaxa) :: lpindexw

#endif

      ! Generate The markers index vector for different structural elements:
      ! In this marker index array the marker points are sorted by structural element.
      !integer ib_MKELEM(ib_nmaxa)

      ! Array of distributed forces:
#if NDIM == 2
      real, save :: ib_Fd(2,ib_nmaxa), ib_Fdaux(2,ib_nmaxa) ! Fstang, Pressure
#elif NDIM == 3
      real, save :: ib_Fd(4,ib_nmaxa), ib_Fdaux(4,ib_nmaxa) ! Fxtang, Fytang, Fztang, Pressure
#endif

      ! Array of aerodynamic elements (segments in 2D or triangles in 3D):
      integer ib_AELEM(ib_nodelema,ib_emaxa)
        
      ! Aerodinamic body type:
      Type ib_abod
         integer :: model, memflag, rigflexflag, lb, mb, elb, emb
         real    :: ron(NDIM)
      end Type ib_abod

      Type(ib_abod), allocatable, dimension(:) :: ib_ABODY

      ! Strings with the names of the bodies files
      ! structural parts
      Type ib_astrtype
        character*60 :: file
      End Type ib_astrtype

      Type(ib_astrtype), allocatable, dimension(:) :: ib_STRA

      ! Aerodynamic elements areas:
      real, save :: ib_AAREAELEM(ib_emaxa)

!----- Variables related to the moving least squares, interpolation functions.
      
!     Interpolation degree for MLS interpolation scheme
!     1 = Linear interpolation and 2 = Quadratic.
      integer, parameter :: ib_interp = 1
      integer, parameter :: ib_npol = NDIM + 1
      integer, parameter :: ib_nderiv = NDIM + 1

!     Number of Eulerian points on the interpolation stencil (5 in 2D, 7 in 3D):      
      integer, parameter :: ib_stencil= 2*NDIM + 1

      integer, save :: ib_ibd

#if NDIM == 2
       real, save :: ib_dutdn(ib_nmaxa),ib_dwtdn(ib_nmaxa),ib_pn(ib_nmaxa), &
                     ib_pb(ib_nmaxa),ib_omb(ib_nmaxa),ib_sb(ib_nmaxa)
#endif

    
!     Arrays for use in the MLS interpolation Scheme:
      real, save :: ib_dsxu(ib_nmaxa),ib_dsyu(ib_nmaxa)
      real, save :: ib_dsxv(ib_nmaxa),ib_dsyv(ib_nmaxa)
      real, save :: ib_dsxu2(ib_nmaxa),ib_dsyv2(ib_nmaxa) 
      integer, save :: ib_ielemu(ib_stencil,ib_nmaxa),ib_jelemu(ib_stencil,ib_nmaxa)
      integer, save :: ib_ielemv(ib_stencil,ib_nmaxa),ib_jelemv(ib_stencil,ib_nmaxa)
      real, save :: ib_UL(ib_nmaxa),ib_VL(ib_nmaxa)
      real, save :: ib_FuL(ib_nmaxa),ib_FvL(ib_nmaxa)
      real, save :: ib_FuL2(ib_nmaxa),ib_indFuL(ib_nmaxa),ib_FvL2(ib_nmaxa),ib_indFvL(ib_nmaxa)
      real, save :: ib_phileu(ib_stencil,ib_nmaxa),ib_philev(ib_stencil,ib_nmaxa)
      real, allocatable, dimension(:) :: ib_xminb,ib_xmaxb,ib_yminb,ib_ymaxb

      

      integer, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS) :: flaguo,flagui
      integer, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,MAXBLOCKS) :: flagvo,flagvi
      integer, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS) :: flagpo,flagpi

      !integer, dimension(:,:,:,:), allocatable :: flaguo,flagui
      !integer, dimension(:,:,:,:), allocatable :: flagvo,flagvi
      !integer, dimension(:,:,:,:), allocatable :: flagpo,flagpi

      real, parameter :: ib_alphax = 1.2
      real, parameter :: ib_alphay = 1.2

      !integer :: ib_bflags(MDIM,MAXBLOCKS)

!     Extended points positions:
      real, save :: ib_xbe(ib_nmaxa),ib_ybe(ib_nmaxa)

      logical,allocatable,dimension(:,:,:,:) :: logvarx,logvary

#if NDIM == 3
      real, save :: ib_dszu(ib_nmaxa)
      real, save :: ib_dszv(ib_nmaxa)
      real, save :: ib_dsxw(ib_nmaxa),ib_dsyw(ib_nmaxa),ib_dszw(ib_nmaxa)
      integer, save :: ib_kelemu(ib_stencil,ib_nmaxa)
      integer, save :: ib_kelemv(ib_stencil,ib_nmaxa)
      integer, save :: ib_kelemw(ib_stencil,ib_nmaxa)
      real, save :: ib_WL(ib_nmaxa)
      real, save :: ib_FwL(ib_nmaxa)
      real, save :: ib_FwL2(ib_nmaxa),ib_indFwL(ib_nmaxa)
      real, save :: ib_philew(ib_stencil,ib_nmaxa)
      real, allocatable, dimension(:) :: ib_zminb,ib_zmaxb

      integer, dimension(:,:,:,:), allocatable :: flagwo,flagwi

      real, parameter :: ib_alphaz = 1.2

!     Extended points positions:
      real, save :: ib_zbe(ib_nmaxa)

      logical,allocatable,dimension(:,:,:,:) :: logvarz      
#endif


      
     real, save :: freq_nat,freq_t,omg,Ro,ao
     real, save,allocatable,dimension(:) :: xo,yo

      integer,save :: ib_meshMe,   ib_meshComm !!, ib_meshNumProcs
      integer,save :: ib_globalMe !!, ib_globalComm, ib_globalNumProcs

end module ImBound_data

