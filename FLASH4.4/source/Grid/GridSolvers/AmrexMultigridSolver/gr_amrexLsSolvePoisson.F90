!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsSolvePoisson
!!
!!  NAME 
!!
!! gr_amrexLsSolvePoisson
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsSolvePoisson()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Poisson equation from the 
!! Amrex Linear Solvers. This is called from  Grid_solvePoisson
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

subroutine gr_amrexLsSolvePoisson ()
  
!   use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
!   use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
!                                Grid_getBlkPtr, Grid_releaseBlkPtr,        &
!                                Grid_getBlkIndexLimits, Grid_getBlkData,   &
!                                Grid_getBlkRefineLevel

  use gr_amrexLsData, ONLY : gr_amrexLs_geom, gr_amrexLs_ba, gr_amrexLs_dm, &
!                                             gr_amrexLs_ascalar, gr_amrexLs_bscalar, &
                                            gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                            gr_amrexLs_max_level, gr_amrexLs_ref_ratio, &
                                            gr_amrexLs_rhs, gr_amrexLs_solution, &
                                            gr_amrexLs_verbose, gr_amrexLs_cg_verbose,        &
                                            gr_amrexLs_max_iter, gr_amrexLs_max_fmg_iter, &
                                            gr_amrexLs_composite_solve, gr_amrexLs_linop_maxorder
!   
  use amrex_fort_module,     ONLY : amrex_real
  use amrex_box_module,     ONLY : amrex_box
!   use amrex_boxarray_module,     ONLY : amrex_boxarray
!   use amrex_distromap_module,     ONLY : amrex_distromap
!   use amrex_geometry_module, ONLY : amrex_geometry, amrex_problo, amrex_probhi
!   use amrex_multifab_module, ONLY : amrex_multifab, amrex_mfiter, &
!                                                 amrex_mfiter_build, amrex_mfiter_destroy
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy
  use amrex_poisson_module, ONLY : amrex_poisson, amrex_poisson_build, amrex_poisson_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic


  
  implicit none
  
#include "Flash.h"
#include "constants.h"   

    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid
    integer :: ilev
    real(amrex_real) :: err
    
  call Timers_start("gr_amrexLsSolvePoisson")     
  print*, "Inside new solve gr_amrexLsSolvePoisson. composite solve=", gr_amrexLs_composite_solve
    if (gr_amrexLs_composite_solve) then

       call amrex_poisson_build(poisson, gr_amrexLs_geom, gr_amrexLs_ba, gr_amrexLs_dm, &
            metric_term=.false., agglomeration=gr_amrexLs_agglomeration, consolidation=gr_amrexLs_consolidation)
       
       call poisson % set_maxorder(gr_amrexLs_linop_maxorder)

!        ! This is a 3d problem with Dirichlet BC
!        call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet], &
!             &                       [amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet])
          call poisson % set_domain_bc([amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic], &
               &                       [amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic])
!	Now solving 2d problem since Flash is not currently set up for 3d
!       call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet], &
!            &                       [amrex_lo_dirichlet, amrex_lo_dirichlet])

       do ilev = 0, gr_amrexLs_max_level
          ! solution multifab's ghost cells at physical boundaries have been set to bc values.
          call poisson % set_level_bc(ilev, gr_amrexLs_solution(ilev))
       end do

       call amrex_multigrid_build(multigrid, poisson)
       call multigrid % set_verbose(gr_amrexLs_verbose)
       call multigrid % set_cg_verbose(gr_amrexLs_cg_verbose)
       call multigrid % set_max_iter(gr_amrexLs_max_iter)
       call multigrid % set_max_fmg_iter(gr_amrexLs_max_fmg_iter)

       print*, "calling multigrid solve, maxlev", gr_amrexLs_max_level
!        print*, "calling multigrid solve, maxlev" gr_amrexLs_max_level
       err = multigrid % solve(gr_amrexLs_solution, gr_amrexLs_rhs, 1.e-10_amrex_real, 0.0_amrex_real)
        print*, err
       call amrex_poisson_destroy(poisson)
       call amrex_multigrid_destroy(multigrid)

    else
       do ilev = 0, gr_amrexLs_max_level

          call amrex_poisson_build(poisson, [gr_amrexLs_geom(ilev)], [gr_amrexLs_ba(ilev)], [gr_amrexLs_dm(ilev)], &
               metric_term=.false., agglomeration=gr_amrexLs_agglomeration, consolidation=gr_amrexLs_consolidation)
          call poisson % set_maxorder(gr_amrexLs_linop_maxorder)

          ! The order of the following set bc calls matters.

          ! This is a 3d problem with Dirichlet BC amrex_lo_periodic
!           call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet], &
!                &                       [amrex_lo_dirichlet, amrex_lo_dirichlet, amrex_lo_dirichlet])
          call poisson % set_domain_bc([amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic], &
               &                       [amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic])
!	Now solving 2d problem since Flash is not currently set up for 3d
!          call poisson % set_domain_bc([amrex_lo_dirichlet, amrex_lo_dirichlet], &
!               &                       [amrex_lo_dirichlet, amrex_lo_dirichlet])
               
          if (ilev > 0) then
             ! use coarse level data to set up bc at corase/fine boundary
             call poisson % set_coarse_fine_bc(gr_amrexLs_solution(ilev-1), gr_amrexLs_ref_ratio)
          end if

          ! Note that to the linear solver, the level is ZERO.  In
          ! this test problem, when lev > 0, solution(lev) is going to
          ! be ignored because fine level grids are completed
          ! surrounded by coarse level.  If fine level grids do touch
          ! phyical domain, the multifab must have bc values at
          ! physical boundaries stored in ghost cells.
          call poisson % set_level_bc(0, gr_amrexLs_solution(ilev))

          call amrex_multigrid_build(multigrid, poisson);
          call multigrid % set_verbose(gr_amrexLs_verbose)
          call multigrid % set_cg_verbose(gr_amrexLs_cg_verbose)
          call multigrid % set_max_iter(gr_amrexLs_max_iter)
          call multigrid % set_max_fmg_iter(gr_amrexLs_max_fmg_iter)

          err = multigrid % solve([gr_amrexLs_solution(ilev)], [gr_amrexLs_rhs(ilev)], 1.e-10_amrex_real, 0.0_amrex_real)

          call amrex_poisson_destroy(poisson)
          call amrex_multigrid_destroy(multigrid)

       end do
    end if
    
  call Timers_stop("gr_amrexLsSolvePoisson")    

  end subroutine gr_amrexLsSolvePoisson