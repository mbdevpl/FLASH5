!!****if* source/Grid/GridMain/AMR/Amrex/gr_amrexInit
!!
!! NAME
!!  gr_amrexInit
!!
!! SYNOPSIS
!!  call gr_amrexInit()
!! 
!! DESCRIPTION
!!  Configure AMReX using FLASH-supplied runtime parameters and for operating in
!!  octree mode.
!!
!!***

subroutine gr_amrexInit()
  use iso_c_binding
  
  use amrex_init_module,           ONLY : amrex_init
  use amrex_amrcore_module,        ONLY : amrex_amrcore_init, &
                                          amrex_init_virtual_functions, &
                                          amrex_max_level, &
                                          amrex_ref_ratio
  use amrex_octree_module,         ONLY : amrex_octree_init
  use amrex_parmparse_module,      ONLY : amrex_parmparse, &
                                          amrex_parmparse_build, &
                                          amrex_parmparse_destroy
  use amrex_geometry_module,       ONLY : amrex_pmask

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use Driver_interface,            ONLY : Driver_abortFlash
  use Grid_data,                   ONLY : gr_geometry, &
                                          gr_domainBC, &
                                          gr_meshMe
  use gr_amrexInterface,           ONLY : gr_initNewLevelCallback, &
                                          gr_makeFineLevelFromCoarseCallback, &
                                          gr_remakeLevelCallback, &
                                          gr_clearLevelCallback, &
                                          gr_markRefineDerefineCallback

  implicit none

#include "constants.h"

  ! From AMReX_CoordSys.H
  ! DEVNOTE: These should be define in AMReX interface.  Name must 
  ! note overlap those used in FLASH
  integer,  parameter :: AMREX_CARTESIAN   = 0
  integer,  parameter :: AMREX_CYLINDRICAL = 1
  integer,  parameter :: AMREX_SPHERICAL   = 2

  type(amrex_parmparse) :: pp_geom
  type(amrex_parmparse) :: pp_amr

  character(len=MAX_STRING_LENGTH) :: buffer = ""

  integer :: verbosity = 0
  integer :: nblockX = 1
  integer :: nblockY = 1
  integer :: nblockZ = 1
  integer :: max_refine = 0
  integer :: nrefs = 1
  integer :: coord_sys = -1
  real    :: xmin = 0.0d0
  real    :: xmax = 1.0d0
  real    :: ymin = 0.0d0
  real    :: ymax = 1.0d0
  real    :: zmin = 0.0d0
  real    :: zmax = 1.0d0
  integer :: is_periodic(MDIM) = 0
  integer :: is_periodic_am(MDIM) = 0

  if(gr_meshMe==MASTER_PE) write(*,*) "[gr_amrexInit] Starting"
 
  !!!!!----- INITIALIZE AMReX & CONFIGURE MANUALLY
  ! Do not parse command line or any file for configuration
 
  ! DEVNOTE: Let FLASH driver construct communicators and then configure AMReX
  ! to use it here.
  call amrex_init(arg_parmparse=.FALSE.)
 
  call amrex_parmparse_build(pp_geom, "geometry")

  select case (gr_geometry)
  case(CARTESIAN)
    coord_sys = AMREX_CARTESIAN
  case(POLAR)
    ! DEVNOTE: FIXME
    call Driver_abortFlash("Not certain how Flash/polar maps only AMReX")
  case(CYLINDRICAL)
    coord_sys = AMREX_CYLINDRICAL
  case(SPHERICAL)
    coord_sys = AMREX_SPHERICAL
  case default
    write(buffer,'(I5)') gr_geometry
    call Driver_abortFlash("Unknown coordinate system type - " &
                           // TRIM(ADJUSTL(buffer)))
  end select
  call pp_geom%add   ("coord_sys", coord_sys)
 
  call RuntimeParameters_get('xmin', xmin)
  call RuntimeParameters_get('xmax', xmax)
  call RuntimeParameters_get('ymin', ymin)
  call RuntimeParameters_get('ymax', ymax)
  call RuntimeParameters_get('zmin', zmin)
  call RuntimeParameters_get('zmax', zmax)
  call pp_geom%addarr("prob_lo", [xmin, ymin, zmin])
  call pp_geom%addarr("prob_hi", [xmax, ymax, zmax])
 
  is_periodic = 0
  if (     (gr_domainBC(LOW,  IAXIS) == PERIODIC) &
      .OR. (gr_domainBC(HIGH, IAXIS) == PERIODIC)) then
    is_periodic(IAXIS) = 1
  end if
  if (     (gr_domainBC(LOW,  JAXIS) == PERIODIC) &
      .OR. (gr_domainBC(HIGH, JAXIS) == PERIODIC)) then
    is_periodic(JAXIS) = 1
  end if
  if (     (gr_domainBC(LOW,  KAXIS) == PERIODIC) &
      .OR. (gr_domainBC(HIGH, KAXIS) == PERIODIC)) then
    is_periodic(KAXIS) = 1
  end if
  call pp_geom%addarr("is_periodic", is_periodic)

  call amrex_parmparse_build(pp_amr, "amr")
  call RuntimeParameters_get("gr_amrex_verbosity", verbosity)
  call pp_amr%add   ("v", verbosity)

  call RuntimeParameters_get("nrefs", nrefs)
  call pp_amr%add   ("regrid_int", nrefs)

  ! AMReX uses 0-based level index set
  call RuntimeParameters_get('lrefine_max', max_refine)
  call pp_amr%add   ("max_level", max_refine - 1)

  call RuntimeParameters_get("nblockx", nBlockX)
  call RuntimeParameters_get("nblocky", nBlockY)
  call RuntimeParameters_get("nblockz", nblockZ)
  call pp_amr%addarr("n_cell", [NXB * nBlockX, &
                                NYB * nBlockY, &
                                NZB * nBlockZ])

  ! Setup AMReX in octree mode
  call pp_amr%add   ("max_grid_size_x",     NXB) 
  call pp_amr%add   ("max_grid_size_y",     NYB) 
  call pp_amr%add   ("max_grid_size_z",     NZB) 
  call pp_amr%add   ("blocking_factor_x", 2*NXB)
  call pp_amr%add   ("blocking_factor_y", 2*NYB)
  call pp_amr%add   ("blocking_factor_z", 2*NZB)
  call pp_amr%add   ("refine_grid_layout", 0)
 
  ! According to Weiqun n_proper=1 is an appropriate setting that will result in
  ! correct nesting.
  call pp_amr%add("n_proper", 1)
  ! We account for data in the guardcells when determining error metric for
  ! entire block.  Therefore, we do not need to tag guardcells.
  call pp_amr%add("n_error_buf", 0)

  ! desctructors not valid for all compilers
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
  call amrex_parmparse_destroy(pp_geom)
  call amrex_parmparse_destroy(pp_amr)
#endif

  call amrex_octree_init()
  call amrex_amrcore_init()
 
  !!!!!----- REGISTER REFINE CALLBACKS WITH AMReX
  call amrex_init_virtual_functions(gr_initNewLevelCallback, &
                                    gr_makeFineLevelFromCoarseCallback, &
                                    gr_remakeLevelCallback, &
                                    gr_clearLevelCallback, &
                                    gr_markRefineDerefineCallback)

  !!!!!----- CONFIRM CORRECT AMReX CONFIGURATION
  if (SIZE(amrex_ref_ratio) /= (max_refine - 1)) then
    call Driver_abortFlash("[gr_amrexInit] AMReX ref_ratio has wrong size") 
  else if (.NOT. ALL(amrex_ref_ratio == 2)) then
    call Driver_abortFlash("[gr_amrexInit] AMReX ref_ratio has wrong values") 
  end if

  ! Check AMReX-controlled BC information
  is_periodic_am(:) = 0
  where (amrex_pmask)  is_periodic_am = 1

  if (is_periodic_am(IAXIS) /= is_periodic(IAXIS)) then
    call Driver_abortFlash("[gr_amrexInit] AMReX does not have correct periodicity in X") 
  end if
#if NDIM >= 2
  if (is_periodic_am(JAXIS) /= is_periodic(JAXIS)) then
    call Driver_abortFlash("[gr_amrexInit] AMReX does not have correct periodicity in Y") 
  end if
#if NDIM == 3
  if (is_periodic_am(KAXIS) /= is_periodic(KAXIS)) then
    call Driver_abortFlash("[gr_amrexInit] AMReX does not have correct periodicity in Z") 
  end if
#endif
#endif
  
  if(gr_meshMe==MASTER_PE) write(*,*) "[gr_amrexInit] Finished"
end subroutine gr_amrexInit

