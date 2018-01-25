!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Driver routine for poisson solvers in the Grid
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!           !DEV: requested bcTypes ignored, always uses GRID_PDE_BND_NEUMANN
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!  
!!  Most implementations only support a limited subset of boundary condition types.
!!  Some implementations require that all significant elements of bcTypes are the same.
!!  (That is always the case for GRID_PDE_BND_ISOLATED.)
!!  Even if an implementation supports combinations of different boundary conditions
!!  on different sides of the domain, the types at left and right sides for the same
!!  axis direction will usually have to be the samme.
!!
!!  Support in some implementations provided with FLASH4:
!!
!!   GridSolvers/Multipole:                  GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/Pfft:                       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.), 
!!                                            GRID_PDE_BND_NEUMANN(hom.)
!!                                            in various combinations,
!!                                            depending on GridSolvers/Pfft subdirectory 
!!                                            (i.e. implementation configured in)
!!                                           (requires UG in pencil shape or Paramesh as Grid)
!!   GridSolvers/Multigrid hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations,
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/HYPRE:                      currently GRID_PDE_BND_NEUMANN only
!!   
!!***

!!-----Do not reorder----!!REORDER(4): solnVec

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
!   use gr_interface,     ONLY : gr_hypreCreateMatrix, gr_hypreComputeB,    &
!                                gr_hypreGridStatus
!   use gr_hypreLocalInterface, ONLY: gr_hypreExchangeFacB
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,        &
                               Grid_getBlkIndexLimits, Grid_getBlkData,   &
                               Grid_getBlkRefineLevel

  use gr_amrexLsData, ONLY : gr_amrexLs_geom, gr_amrexLs_ba, gr_amrexLs_dm, &
                                            gr_amrexLs_ascalar, gr_amrexLs_bscalar, &
                                            gr_amrexLs_max_level, gr_amrexLs_prob_type, &
                                            gr_amrexLs_exact_solution, gr_amrexLs_rhs, gr_amrexLs_solution
!   use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
!                                gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
!                                gr_hypreUseFloor
!   
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET

use mytest_module, ONLY : init, finalize, solve, write_plotfile
use gr_amrexLsInterface, ONLY : gr_amrexLsInitPoisson, gr_amrexLsInitGeom, &
                                          gr_amrexLsInitGeom

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW  
  
  integer :: mylevel, mypart, var, ii,i,j,k, ierr
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID, lb
  real, allocatable :: cellVolumes(:,:,:)
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  logical :: mask(NUNK_VARS), savedUseFloor

  integer :: iFactorB
  integer :: iFactorA
  real    :: dt, theta

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  integer :: localbcTypes(6) 
!!$    character(len=32) :: matfile

  
  call Timers_start("Grid_solvePoisson")    
  print*, "Inside Grid_solvePoisson using amrex solvers"
  !!*********************----------------*******************!!
  !allocate data
    allocate(gr_amrexLs_geom(0:gr_amrexLs_max_level))
    allocate(gr_amrexLs_ba(0:gr_amrexLs_max_level))
    allocate(gr_amrexLs_dm(0:gr_amrexLs_max_level))
    allocate(gr_amrexLs_solution(0:gr_amrexLs_max_level))
    allocate(gr_amrexLs_rhs(0:gr_amrexLs_max_level))
    allocate(gr_amrexLs_exact_solution(0:gr_amrexLs_max_level))
    if (gr_amrexLs_prob_type .eq. 2) then
       call Driver_abortFlash('Abec problem not implemnted yet!! Program will abort!!!!')
    end if
    call gr_amrexLsInitGeom ()
    call gr_amrexLsInitGrid()
    call gr_amrexLsInitMf()
    if (gr_amrexLs_prob_type .eq. 1) then
       call gr_amrexLsInitPoisson(gr_amrexLs_geom,gr_amrexLs_solution, gr_amrexLs_rhs, gr_amrexLs_exact_solution)
!        call gr_amrexLsSolvePoisson()
    else
!        call init_prob_abeclaplacian(gr_amrexLs_geom,solution, gr_amrexLs_rhs, gr_amrexLs_exact_solution, &
!        gr_amrexLs_acoef, gr_amrexLs_bcoef, gr_amrexLs_ascalar, gr_amrexLs_bscalar)
!          TODO :: abec solve
         call Driver_abortFlash('Abec problem not implemnted yet!! Program will abort!!!!')
    end if
!     call gr_amrexLsInitPoisson (gr_amrexLs_geom, gr_amrexLs_solution, gr_amrexLs_rhs, gr_amrexLs_exact_solution)
    
    
    
  ! allocate data
  call init()
  call solve()
  call finalize()


!   call amrex_finalize()
  
  
  
  !!*********************----------------*******************!!
  
  
  
  call Timers_stop("Grid_solvePoisson") 
  
  return
end subroutine Grid_solvePoisson
