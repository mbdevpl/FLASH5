!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsData
!!
!! NAME
!!
!!  gr_AmrexLsData
!!
!! SYNOPSIS
!!  use gr_AmrexLsData
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for the HYPRE implementation
!!
!!  
!!***

#include "Flash.h"

module gr_amrexLsData
  
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
  use amrex_multifab_module, ONLY : amrex_multifab
  use amrex_boxarray_module,     ONLY : amrex_boxarray
  use amrex_distromap_module,     ONLY : amrex_distromap
  use amrex_geometry_module, ONLY : amrex_geometry
  use amrex_fort_module,     ONLY : amrex_real
  
  implicit none

  ! parameters

  integer, save :: gr_amrexLs_max_level
  integer, save :: gr_amrexLs_ref_ratio = 2
  integer, save :: gr_amrexLs_n_cell
  integer, save :: gr_amrexLs_max_grid_size

  logical, save :: gr_amrexLs_composite_solve

  ! prob_type 1 here is Poisson with homogeneous Dirichlet boundary.
  ! prob_type 2 here is ABecLaplacian with homogeneous Neumann boundary.
  integer, save :: gr_amrexLs_prob_type = 1

  integer, save :: gr_amrexLs_verbose = 2
  integer, save :: gr_amrexLs_cg_verbose = 0
  integer, save :: gr_amrexLs_max_iter = 100
  integer, save :: gr_amrexLs_max_fmg_iter = 0
  integer, save :: gr_amrexLs_linop_maxorder = 2
  logical, save :: gr_amrexLs_agglomeration = .true.
  logical, save :: gr_amrexLs_consolidation = .true.

  ! data
  type(amrex_geometry), allocatable, save :: gr_amrexLs_geom(:)
  type(amrex_boxarray), allocatable, save :: gr_amrexLs_ba(:)
  type(amrex_distromap), allocatable, save :: gr_amrexLs_dm(:)

  type(amrex_multifab), allocatable, save :: gr_amrexLs_solution(:)
  type(amrex_multifab), allocatable, save :: gr_amrexLs_rhs(:)
  type(amrex_multifab), allocatable, save :: gr_amrexLs_exact_solution(:)
  type(amrex_multifab), allocatable, save :: gr_amrexLs_acoef(:)
  type(amrex_multifab), allocatable, save :: gr_amrexLs_bcoef(:)
  real(amrex_real), save :: gr_amrexLs_ascalar, gr_amrexLs_bscalar
  
end module gr_amrexLsData
