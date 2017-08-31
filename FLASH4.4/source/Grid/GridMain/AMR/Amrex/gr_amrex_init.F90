subroutine gr_amrex_init()
  use iso_c_binding
  
  use amrex_fort_module,      ONLY : wp => amrex_real
  use amrex_amr_module,       ONLY : amrex_init, &
                                     amrex_amrcore_init, &
                                     amrex_init_virtual_functions, &
                                     amrex_init_from_scratch, &
                                     amrex_max_level
  use amrex_parmparse_module, ONLY : amrex_parmparse, &
                                     amrex_parmparse_build, &
                                     amrex_parmparse_destroy
  use amrex_octree_module,    ONLY : amrex_octree_init
  use amrex_interfaces,       ONLY : gr_makeNewLevelFromScratch, &
                                     gr_makeNewLevelFromCoarse, &
                                     gr_remakeLevel, &
                                     gr_clearLevel, &
                                     gr_markRefineDerefine
  use physicaldata,           ONLY : unk

  use Driver_interface,       ONLY : Driver_abortFlash
  use Grid_data,              ONLY : gr_verbosity, &
                                     gr_geometry, &
                                     gr_imin, gr_imax, &
                                     gr_jmin, gr_jmax, &
                                     gr_kmin, gr_kmax, &
                                     gr_nBlockX, gr_nBlockY, gr_nBlockZ, &
                                     gr_maxRefine, &
                                     gr_nrefs

  implicit none

#include "constants.h"

  ! From AMReX_CoordSys.H
  ! DEVNOTE: These should be define in AMReX interface.  Name must 
  ! note overlap those used in FLASH
  integer,  parameter :: AMREX_CARTESIAN   = 0
  integer,  parameter :: AMREX_CYLINDRICAL = 1
  integer,  parameter :: AMREX_SPHERICAL   = 2

  real(wp), parameter :: T_INIT = 0.0_wp

  type(amrex_parmparse) :: pp_geom
  type(amrex_parmparse) :: pp_amr

  integer       :: coord_sys = -1
  character(80) :: buffer = ""

  write(*,*) "[gr_amrex_init] Starting"
 
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
  call pp_geom%addarr("prob_lo", [gr_imin, gr_jmin, gr_kmin])
  call pp_geom%addarr("prob_hi", [gr_imax, gr_jmax, gr_kmax])
  call pp_geom%addarr("is_periodic", [1, 1, 1])

  ! DEVNOTE: max_grid must be a multiple of blocking_factor.  Error checking in
  ! AMReX or here?
  call amrex_parmparse_build(pp_amr, "amr")
  call pp_amr%add   ("v", gr_verbosity)
  call pp_amr%add   ("max_level", gr_maxRefine)
  call pp_amr%add   ("regrid_int", gr_nrefs)
  call pp_amr%addarr("n_cell", [NXB * gr_nBlockX, &
                                NYB * gr_nBlockY, &
                                NZB * gr_nBlockZ])
  call pp_amr%add   ("max_grid_size", 4)
  call pp_amr%add   ("blocking_factor", 2)
!  call pp_amr%add   ("n_proper", )
!  call pp_amr%add   ("grid_eff", ._wp)
!  call pp_amr%add   ("n_error_buf", )
  call pp_amr%add   ("refine_grid_layout", 0)

    ! Finalization not supported by all compilers.  Look at AMReX iterators for 
    ! preprocessor directives to isolate these.
!  call amrex_parmpase_destroy(pp_geom)
!  call amrex_parmpase_destroy(pp_amr)

  call amrex_octree_init()
  call amrex_amrcore_init()
 
  !!!!!----- REGISTER REFINE CALLBACKS WITH AMReX
  call amrex_init_virtual_functions(gr_makeNewLevelFromScratch, &
                                    gr_makeNewLevelFromCoarse, &
                                    gr_remakeLevel, &
                                    gr_clearLevel, &
                                    gr_markRefineDerefine)

  !!!!!----- ALLOCATE DATA STRUCTURES
  ! multifabs
  allocate(unk(0:amrex_max_level))

  ! Setup grids and initialize the data
  call amrex_init_from_scratch(T_INIT)
end subroutine gr_amrex_init

