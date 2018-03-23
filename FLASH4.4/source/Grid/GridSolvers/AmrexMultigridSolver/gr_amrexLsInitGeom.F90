!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsInitGeom
!!
!!  NAME 
!!
!! gr_amrexLsInitGeom
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsInitGeom()
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

subroutine gr_amrexLsInitGeom ()
  
!   use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_amrexLsData, ONLY : gr_amrexLs_geom, gr_amrexLs_n_cell, gr_amrexLs_ref_ratio, &
!                                             gr_amrexLs_ascalar, gr_amrexLs_bscalar, &
!                                             gr_amrexLs_max_level, gr_amrexLs_prob_type, &
                                            gr_amrexLs_max_level
  use amrex_fort_module,     ONLY : amrex_real
  use amrex_box_module,     ONLY : amrex_box
  use amrex_geometry_module,    ONLY : amrex_geometry_set_coord_sys, amrex_geometry_set_prob_domain, &
                                                        amrex_geometry_set_periodic, amrex_geometry_build
  implicit none
  
#include "Flash.h"
#include "constants.h"   
  
    integer :: ilev
    type(amrex_box) :: domain

    call Timers_start("gr_amrexLsInitGeom")     
  
    call amrex_geometry_set_coord_sys(0)  ! Cartesian
    call amrex_geometry_set_prob_domain([0._amrex_real,0._amrex_real,0._amrex_real], &
         &                              [1._amrex_real,1._amrex_real,1._amrex_real])
    call amrex_geometry_set_periodic([.false., .false., .false.])
    domain = amrex_box([0,0,0], [gr_amrexLs_n_cell-1,gr_amrexLs_n_cell-1,gr_amrexLs_n_cell-1])
    do ilev = 0, gr_amrexLs_max_level
       call amrex_geometry_build(gr_amrexLs_geom(ilev), domain)
       call domain % refine(gr_amrexLs_ref_ratio)
    end do  
  
  call Timers_stop("gr_amrexLsInitGeom")    

end subroutine gr_amrexLsInitGeom
