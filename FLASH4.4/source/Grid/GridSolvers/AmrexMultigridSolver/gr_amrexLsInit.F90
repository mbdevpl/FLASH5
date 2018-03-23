!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsInit
!!
!! NAME
!!
!!  gr_amrexLsInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_amrexLsInit()
!!
!! Description
!!
!!  Initializes local data for Unit AmrexMultigridSolver defined in Module gr_amrexLsData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  AmrexMultigridSolver.
!!
!! ARGUMENTS
!!
!!  none  
!!***

subroutine gr_amrexLsInit()
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use gr_amrexLsData, ONLY : gr_amrexLs_max_level, gr_amrexLs_ref_ratio,           &
                           gr_amrexLs_n_cell,                                      & 
                           gr_amrexLs_max_grid_size, gr_amrexLs_composite_solve,    & 
                           gr_amrexLs_prob_type, gr_amrexLs_cg_verbose,        &
                           gr_amrexLs_max_iter, gr_amrexLs_max_fmg_iter, &
                           gr_amrexLs_linop_maxorder, gr_amrexLs_verbose,      &
                           gr_amrexLs_agglomeration, gr_amrexLs_consolidation
  
  implicit none 

#include "constants.h"
print*, "reading params in gr_amrexLsInit"
  call RuntimeParameters_get("gr_amrexLs_max_level",   gr_amrexLs_max_level)
  call RuntimeParameters_get("gr_amrexLs_ref_ratio", gr_amrexLs_ref_ratio)
  call RuntimeParameters_get("gr_amrexLs_n_cell", gr_amrexLs_n_cell)
  call RuntimeParameters_get("gr_amrexLs_max_grid_size", gr_amrexLs_max_grid_size)
  
  call RuntimeParameters_get("gr_amrexLs_prob_type", gr_amrexLs_prob_type)
  
  call RuntimeParameters_get("gr_amrexLs_verbose",    gr_amrexLs_verbose)
  call RuntimeParameters_get("gr_amrexLs_cg_verbose",    gr_amrexLs_cg_verbose)

  call RuntimeParameters_get("gr_amrexLs_max_iter",    gr_amrexLs_max_iter)
  call RuntimeParameters_get("gr_amrexLs_max_fmg_iter",    gr_amrexLs_max_fmg_iter)
  call RuntimeParameters_get("gr_amrexLs_linop_maxorder",  gr_amrexLs_linop_maxorder)
  call RuntimeParameters_get("gr_amrexLs_agglomeration", gr_amrexLs_agglomeration)
  call RuntimeParameters_get("gr_amrexLs_consolidation",   gr_amrexLs_consolidation)
  
end subroutine gr_amrexLsInit
