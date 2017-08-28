subroutine gr_amrex_init()
  use iso_c_binding
  
  use amrex_fort_module,      ONLY : wp => amrex_real
  use amrex_amr_module,       ONLY : amrex_init, &
                                     amrex_amrcore_init, &
                                     amrex_init_virtual_functions, &
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

!  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"

  type(amrex_parmparse) :: pp
  type(amrex_parmparse) :: pp_geom
  type(amrex_parmparse) :: pp_amr
 
  !!!!!----- INITIALIZE AMReX & CONFIGURE MANUALLY
  ! Do not parse command line or any file for configuration
  
  ! DEVNOTE: Let FLASH driver construct communicators and then configure AMReX
  ! to use it here.
  call amrex_init(arg_parmparse=.FALSE.)
  
  ! DEVNOTE: These values should be taken from the Runtime FLASH facility
  call amrex_parmparse_build(pp)
  call pp%add("max_step", 1000000)
  call pp%add("stop_time", 2.0_wp)

  ! DEVNOTE: Can we use the C-preprocessor AMREX_D_DECL macro here?  Do we want
  ! to use it?
  call amrex_parmparse_build(pp_geom, "geometry")
  call pp_geom%add   ("coord_sys", CARTESIAN)
  call pp_geom%addarr("prob_lo", [0.0_wp, 0.0_wp, 0.0_wp])
  call pp_geom%addarr("prob_hi", [1.0_wp, 1.0_wp, 1.0_wp])
  call pp_geom%addarr("is_periodic", [1, 1, 1])

  ! DEVNOTE: max_grid must be a multiple of blocking_factor.  Error checking in
  ! AMReX or here?
  call amrex_parmparse_build(pp_amr, "amr")
  call pp_amr%add   ("v", 1)
  call pp_amr%add   ("max_level", 2)
  call pp_amr%add   ("regrid_int", 2)
  call pp_amr%addarr("n_cell", [8, 8, 8])
  call pp_amr%add   ("max_grid_size", 4)
  call pp_amr%add   ("blocking_factor", 2)
!  call pp_amr%add   ("n_proper", )
!  call pp_amr%add   ("grid_eff", ._wp)
!  call pp_amr%add   ("n_error_buf", )
  call pp_amr%add   ("refine_grid_layout", 0)

    ! Finalization not supported by all compilers.  Look at AMReX iterators for 
    ! preprocessor directives to isolate these.
!  call amrex_parmpase_destroy(pp)
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
end subroutine gr_amrex_init

