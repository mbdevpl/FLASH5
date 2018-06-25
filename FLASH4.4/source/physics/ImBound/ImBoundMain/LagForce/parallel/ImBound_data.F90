!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ImBound_data
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

      ! MPI Related:
      integer,save :: ib_meshMe,   ib_meshComm !!, ib_meshNumProcs
      integer,save :: ib_globalMe !!, ib_globalComm, ib_globalNumProcs

      ! Timestep:
      real, save :: ib_dt
  
      ! Kinematic viscosity:
      real, save :: ib_nu = 0.

      ! Maximum iteration of IB forcing
      integer, save :: ib_maxIterForcing

      ! Number of bodies
      integer, save :: ib_nbd
      integer, save :: ib_ibd

      ! Maximum number of marker points in the system
      integer, save :: ib_nnoda

      ! Aerodinamic surfaces elements
      integer, save :: ib_nela

      ! Flag to note if the block posseses any Lagrangian Markers,
      ! Use for refinement-derefinement:
      logical, save :: ib_BlockMarker_flag(MAXBLOCKS)

!----- Variables related to the moving least squares, interpolation functions.
      
      ! Interpolation degree for MLS interpolation scheme
      ! 1 = Linear interpolation and 2 = Quadratic.
      integer, parameter :: ib_interp = 1
      integer, parameter :: ib_npol = NDIM + 1
      integer, parameter :: ib_nderiv = NDIM + 1

      ! Number of Eulerian points on the interpolation stencil (5 in 2D, 7 in 3D):      
      integer, parameter :: ib_stencil= 2*NDIM + 1

      ! Scaling factors for domain of influence of Lagrangian Markers
      real, parameter :: ib_alphax  = 1.2
      real, parameter :: ib_alphay  = 1.2
      real, parameter :: ib_alphaz  = 1.2
      real, parameter :: distFactor = 0.6

!!$      integer, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS) :: flaguo,flagui
!!$      integer, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,MAXBLOCKS) :: flagvo,flagvi
!!$      integer, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)   :: flagpo,flagpi

      ! 2D Cylinder runs Specific data - To be Deleted.
      real, save :: freq_nat,freq_t,omg,Ro,ao
      real, save,allocatable,dimension(:) :: xo,yo,zo

end module ImBound_data

