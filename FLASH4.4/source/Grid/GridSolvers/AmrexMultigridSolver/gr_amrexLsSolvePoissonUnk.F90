!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsSolvePoissonUnk
!!
!!  NAME 
!!
!! gr_amrexLsSolvePoissonUnk
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsSolvePoissonUnk()
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
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Solver settings used ?? DESCRIPTION??
!!
!!***

subroutine gr_amrexLsSolvePoissonUnk ()
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy
  use amrex_poisson_module, ONLY : amrex_poisson, amrex_poisson_build, amrex_poisson_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic
  use amrex_amr_module, ONLY : amrex_geom, amrex_get_finest_level
  use amrex_fort_module,     ONLY : amrex_real
  use gr_amrexLsData, ONLY : gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                    gr_amrexLs_linop_maxorder, gr_amrexLs_verbose, gr_amrexLs_cg_verbose, &
                                    gr_amrexLs_max_iter, gr_amrexLs_max_fmg_iter
  use gr_physicalMultifabs,  ONLY : unk

  implicit none
  
    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid
    integer :: ilev, maxLevel
    real(amrex_real) :: err
#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("gr_amrexLsInitPoissonUnk")
!   Build poisson object with the geometry amrex_geom, boxarray unk%ba  and distromap unk%dm
       maxLevel = amrex_get_finest_level() !! TODO :: Check with Jared if this is the best wat to get max level of current 
                                                                   !! grid. There is likely a flash ssubroutine for same Grid_getMaxRefinement?
       call amrex_poisson_build(poisson, amrex_geom(0:maxLevel), unk%ba, unk%dm, &
            metric_term=.false., agglomeration=gr_amrexLs_agglomeration, consolidation=gr_amrexLs_consolidation)
       
       call poisson % set_maxorder(gr_amrexLs_linop_maxorder)

!        ! This is a 3d problem with Periodic BC
          call poisson % set_domain_bc([amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic], &
               &                       [amrex_lo_periodic, amrex_lo_periodic, amrex_lo_periodic])

       print* , "Refinement level max",maxLevel
       do ilev = 0, maxLevel
          ! solution multifab's ghost cells at physical boundaries have been set to bc values.
          call poisson % set_level_bc(ilev, unk(ilev))  !!Not sure if this is correct because unk multifab has other vars
       end do

       call amrex_poisson_destroy(poisson)
       
  call Timers_stop("gr_amrexLsInitPoissonUnk")
end subroutine gr_amrexLsSolvePoissonUnk
