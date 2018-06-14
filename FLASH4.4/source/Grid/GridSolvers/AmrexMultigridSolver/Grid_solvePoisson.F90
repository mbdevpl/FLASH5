!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/Grid_solvePoisson
!!
!!  NAME 
!!
!! Grid_solvePoisson
!!
!!  SYNOPSIS
!!
!!  call Grid_solvePoisson()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Poisson equation from the 
!! Amrex Linear Solvers using the variables from Unk multifab
!! for rhs and unknown phi 
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
!!
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Currently, solver only works for GRID_PDE_BND_PERIODIC i.e. periodic boundary conditions
!!  Other BCs to be implemented later
!!  Relative and absolute tolerences for multigrid sovle - 1.e-10, 0.0
!!
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy
  use amrex_poisson_module, ONLY : amrex_poisson, amrex_poisson_build, amrex_poisson_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic, amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_amr_module, ONLY : amrex_geom, amrex_get_finest_level, amrex_max_level
  use amrex_fort_module,     ONLY : amrex_real
  use gr_amrexLsData, ONLY : gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                    gr_amrexLs_linop_maxorder, gr_amrexLs_verbose, gr_amrexLs_cg_verbose, &
                                    gr_amrexLs_max_iter, gr_amrexLs_max_fmg_iter
  use gr_physicalMultifabs,  ONLY : unk
  !!
  use amrex_multifab_module, ONLY : amrex_multifab, amrex_multifab_destroy, amrex_multifab_build_alias

  implicit none
  
    integer, intent(in)    :: iSoln, iSrc
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(inout)    :: poisfact
    integer                :: amrexPoissonBcTypes(6)
    integer                :: i
    
    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid
    integer :: ilev, maxLevel
    real(amrex_real) :: err
    type(amrex_multifab), allocatable, save :: solution(:)
    type(amrex_multifab), allocatable, save :: rhs(:)

#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("Grid_solvePoisson")
     maxLevel = amrex_get_finest_level() !! TODO :: Check with Jared if this is the best wat to get max level of current 
                                                                   !! grid. There is likely a flash ssubroutine for same Grid_getMaxRefinement?
!   Allocate space for multifab array storing phi (solution) and rhs
    allocate(solution(0:maxLevel))
    allocate(rhs(0:maxLevel))

!     Create alias from multifab unk to rhs and solution with respective components in unk
    do ilev = 0, maxLevel
        call amrex_multifab_build_alias(solution(ilev), unk(ilev), iSoln, 1)
        call amrex_multifab_build_alias(rhs(ilev), unk(ilev), iSrc, 1)
        call solution(ilev)%setVal(0.0_amrex_real)
    end do
!   Build poisson object with the geometry amrex_geom, boxarray unk%ba  and distromap unk%dm
       call amrex_poisson_build(poisson, amrex_geom(0:maxLevel), rhs%ba, rhs%dm, &
            metric_term=.false., agglomeration=gr_amrexLs_agglomeration, consolidation=gr_amrexLs_consolidation)
       
       call poisson % set_maxorder(gr_amrexLs_linop_maxorder)

!  Select BCs to send to AMReX poisson solver
     do i=1,6
       select case (bcTypes(i))
       case (GRID_PDE_BND_PERIODIC)
          amrexPoissonBcTypes(i)=amrex_lo_periodic
       case (GRID_PDE_BND_NEUMANN)
          amrexPoissonBcTypes(i)=amrex_lo_neumann
       case (GRID_PDE_BND_DIRICHLET)
          amrexPoissonBcTypes(i)=amrex_lo_dirichlet
       case default
          call Driver_abortFlash('Only periodic BC implemented for AMReX poisson solver!')
       end select
     end do
     call poisson % set_domain_bc([amrexPoissonBcTypes(1),amrexPoissonBcTypes(3),amrexPoissonBcTypes(5)], &
          &                       [amrexPoissonBcTypes(2),amrexPoissonBcTypes(4),amrexPoissonBcTypes(6)])

       do ilev = 0, maxLevel
! solution multifab's ghost cells at physical boundaries have been set to bc values.
          call poisson % set_level_bc(ilev, solution(ilev))
       end do

       call amrex_multigrid_build(multigrid, poisson)
       call multigrid % set_verbose(gr_amrexLs_verbose)
       call multigrid % set_cg_verbose(gr_amrexLs_cg_verbose)
       call multigrid % set_max_iter(gr_amrexLs_max_iter)
       call multigrid % set_max_fmg_iter(gr_amrexLs_max_fmg_iter)

       print*, "Calling multigrid solve, maxlev", maxLevel
       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)
        print*, err
!      !!Finalize objects
        do ilev = 0, maxLevel
            call amrex_multifab_destroy(solution(ilev))
            call amrex_multifab_destroy(rhs(ilev))
        end do
       call amrex_multigrid_destroy(multigrid)
       call amrex_poisson_destroy(poisson)
       
  call Timers_stop("Grid_solvePoisson")
end subroutine Grid_solvePoisson
