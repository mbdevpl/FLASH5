!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsInitGrid
!!
!!  NAME 
!!
!! gr_amrexLsInitGrid
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsInitGrid()
!!
!!
!!  DESCRIPTION 
!! This routine sets up the Poisson solver from the 
!! Amrex Linear Solvers. 
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Solver settings used ?? DESCRIPTION??
!!
!!***

subroutine gr_amrexLsInitGrid ()
  
!   use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_amrexLsData, ONLY : gr_amrexLs_geom, gr_amrexLs_ba, &
                                            gr_amrexLs_max_level, gr_amrexLs_max_grid_size, &
                                            gr_amrexLs_n_cell, gr_amrexLs_ref_ratio
  use amrex_boxarray_module,     ONLY : amrex_boxarray_build
  use amrex_box_module,     ONLY : amrex_box
!   use amrex_distromap_module,     ONLY : amrex_distromap
!   use amrex_geometry_module, ONLY : amrex_geometry
  implicit none
  
#include "Flash.h"
#include "constants.h"   

    integer :: ilev
    type(amrex_box) :: dom
  
    call Timers_start("gr_amrexLsInitGrid")     
  
    dom = gr_amrexLs_geom(0) % domain
    do ilev = 0, gr_amrexLs_max_level
       call amrex_boxarray_build(gr_amrexLs_ba(ilev), dom)
       call gr_amrexLs_ba(ilev) % maxSize(gr_amrexLs_max_grid_size)
       call dom % grow(-gr_amrexLs_n_cell/4)     ! fine level cover the middle of the coarse domain
       call dom % refine(gr_amrexLs_ref_ratio)
    end do
  
  call Timers_stop("gr_amrexLsInitGrid")    

end subroutine gr_amrexLsInitGrid
