!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsInitPoisson
!!
!!  NAME 
!!
!! gr_amrexLsInitPoisson
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsInitPoisson()
!!
!!
!!  DESCRIPTION 
!! This routine sets up the Poisson solver from the 
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

subroutine gr_amrexLsInitPoisson (gr_amrexLs_geom, gr_amrexLs_solution, gr_amrexLs_rhs, gr_amrexLs_exact_solution)
  
!   use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
!   use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getListOfBlocks, &
!                                Grid_getBlkPtr, Grid_releaseBlkPtr,        &
!                                Grid_getBlkIndexLimits, Grid_getBlkData,   &
!                                Grid_getBlkRefineLevel

!   use gr_amrexLsData, ONLY : gr_amrexLs_geom, gr_amrexLs_ba, gr_amrexLs_dm, &
!                                             gr_amrexLs_ascalar, gr_amrexLs_bscalar, &
!                                             gr_amrexLs_max_level, gr_amrexLs_prob_type, &
!                                             gr_amrexLs_exact_solution, gr_amrexLs_rhs, gr_amrexLs_solution
! !   
  use amrex_fort_module,     ONLY : amrex_real
  use amrex_box_module,     ONLY : amrex_box
!   use amrex_boxarray_module,     ONLY : amrex_boxarray
!   use amrex_distromap_module,     ONLY : amrex_distromap
  use amrex_geometry_module, ONLY : amrex_geometry, amrex_problo, amrex_probhi
  use amrex_multifab_module, ONLY : amrex_multifab, amrex_mfiter, &
                                                amrex_mfiter_build, amrex_mfiter_destroy
 

  
  implicit none
  
#include "Flash.h"
#include "constants.h"   
  
    type(amrex_geometry), intent(in   ) :: gr_amrexLs_geom(0:)
    type(amrex_multifab), intent(inout) :: gr_amrexLs_solution(0:)
    type(amrex_multifab), intent(inout) :: gr_amrexLs_rhs(0:)
    type(amrex_multifab), intent(inout) :: gr_amrexLs_exact_solution(0:)

    integer :: ilev
    integer :: rlo(4), rhi(4), elo(4), ehi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: prhs, pexact

  call Timers_start("gr_amrexLsInitPoisson")     
  
    do ilev = 0, size(gr_amrexLs_rhs)-1
       !$omp parallel private(rlo,rhi,elo,ehi,bx,mfi,prhs,pexact)
       call amrex_mfiter_build(mfi, gr_amrexLs_rhs(ilev), tiling=.true.)
       
       do while (mfi%next())
          bx = mfi%tilebox()
          prhs   =>            gr_amrexLs_rhs(ilev) % dataptr(mfi)
          pexact => gr_amrexLs_exact_solution(ilev) % dataptr(mfi)
          rlo = lbound(prhs)
          rhi = ubound(prhs)
          elo = lbound(pexact)
          ehi = ubound(pexact)
          call actual_init_poisson(bx%lo, bx%hi, prhs, rlo(1:3), rhi(1:3), pexact, elo(1:3), ehi(1:3), &
               amrex_problo, amrex_probhi, gr_amrexLs_geom(ilev)%dx)
       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

       ! This will be used to provide bc and initial guess for the solver.
       call gr_amrexLs_solution(ilev)%setVal(0.0_amrex_real)
    end do
    
  call Timers_stop("gr_amrexLsInitPoisson")    
  
  return

contains
    subroutine actual_init_poisson (lo, hi, rhs, rlo, rhi, exact, elo, ehi, prob_lo, prob_hi, dx)
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, elo, ehi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) :: x, y, z
    real(amrex_real), parameter :: tpi =  8.d0*atan(1.0)
    real(amrex_real), parameter :: fpi = 16.d0*atan(1.0)
    real(amrex_real), parameter :: fac = tpi*tpi*3.d0

    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
             exact(i,j,k) = 1.d0 * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &      + .25d0 * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
                
             rhs(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &       -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
          end do
       end do
    end do

  end subroutine actual_init_poisson
  
end subroutine gr_amrexLsInitPoisson
