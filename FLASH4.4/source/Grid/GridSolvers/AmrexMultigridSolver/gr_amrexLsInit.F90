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
!!  Initializes local data for Unit HYPRE defined in Module gr_hypreData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  HYPRE.
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
! #include "gr_hypre_fortran_message.h"
!  
!   character(len=MAX_STRING_LENGTH) :: solver
!   character(len=MAX_STRING_LENGTH) :: pc
!   character(len=22) :: msgBuf(3,4)
!   character(len=HYPRE_MSG_LEN) :: hypreMsg
!   integer :: i,j
!   
!   real :: ssol
  
!   integer, save :: gr_amrexLs_max_level = 1
!   integer, save :: gr_amrexLs_ref_ratio = 2
!   integer, save :: gr_amrexLs_n_cell = 128
!   integer, save :: gr_amrexLs_max_grid_size = 64
! 
!   logical, save :: gr_amrexLs_composite_solve = .true.
! 
!   ! prob_type 1 here is Poisson with homogeneous Dirichlet boundary.
!   ! prob_type 2 here is ABecLaplacian with homogeneous Neumann boundary.
!   integer, save :: gr_amrexLs_prob_type = 1
! 
!   integer, save :: gr_amrexLs_verbose = 2
!   integer, save :: gr_amrexLs_cg_verbose = 0
!   
!   integer, save :: gr_amrexLs_max_iter = 100
!   integer, save :: gr_amrexLs_max_fmg_iter = 0
!   integer, save :: gr_amrexLs_linop_maxorder = 2
!   logical, save :: gr_amrexLs_agglomeration = .true.
!   logical, save :: gr_amrexLs_consolidation = .true.

!   call PhysicalConstants_get("speed of light",gr_speedlt)

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
  
  ! Write version info about the HYPRE library used to the log file.
!   call gr_hypre_fortran_message(hypreMsg, HYPRE_MSG_LEN)
!   call Logfile_stamp(hypreMsg , '[gr_hypreInit]')
! 
!   gr_hypreGridIsSetUp = .FALSE.   
!   call RuntimeParameters_get("nbegin", gr_hypreNStep)
!   gr_hypreNStep = gr_hypreNStep - 1
! 
!   ! Anisotropic diffusion, probably to represent anisotropic thermal conduction
!   call RuntimeParameters_get("gr_hypreAnisoDiffusion",   gr_hypreAnisoDiffusion)
  
end subroutine gr_amrexLsInit
