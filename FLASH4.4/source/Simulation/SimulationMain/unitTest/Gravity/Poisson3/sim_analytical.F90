!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/sim_analytical
!!
!! NAME
!!
!!  sim_analytical
!!
!! SYNOPSIS
!!
!!  sim_analytical(integer(IN) :: tileDesc)
!!
!! DESCRIPTION
!!
!!  Provides analytical solution for the Maclaurin spheroid problem.
!!  Analytical solution given in Ricker, P. 2006 "A Direct Multigrid Poisson Solver for Oct-tree adaptive meshes"
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  tileDesc -- current grid block
!!
!! NOTES
!!
!!  There is a lovely treatment of this topic in a great web page
!!   http://www.phys.lsu.edu/astro/H_Book.current/Applications/Structure/Nonaxisymmetric/Homo_Ellipsoids/OIBT.shtml
!!
!!***

subroutine sim_analytical(tileDesc)

  use Simulation_data, ONLY:  sim_Newton, sim_pi, &
       &  sim_e, sim_a1, sim_a3, sim_a1inv, sim_a3inv, &
       &  sim_xctr, sim_yctr, sim_zctr, &
       &  sim_density, sim_nsubinv, sim_nsubzones, sim_initGeometry

  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_tile,      ONLY : Grid_tile_t
  
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Passed variables
  type(Grid_tile_t), intent(IN) :: tileDesc

  ! local variables
  integer :: lo(1:MDIM)
  integer :: hi(1:MDIM)

  integer       :: sizeX, sizeY, sizeZ
  integer       :: k, kk, i, ii, j, jj
  real          :: a1squared, a3squared, a1squaredInv, a3squaredInv
  real          :: AA1ch, AA3ch, Ich, distTerm
  real          :: DD, Mass, phiSum
  real, dimension(MDIM)    :: deltas
  real          :: dx, dy, dz, xx, yy, zz, dzz, dyy, dxx
  real          :: xdist, ydist, zdist, rInv, rFunc, r2, z2
  real          :: lambda, h, hinv, H0, atan_h, h1, h3

  real, dimension(:), allocatable :: xLeft, yLeft, zLeft

  real, pointer :: solnData(:,:,:,:)

  nullify(solnData)

  ! Convenience variables
  a1squared = sim_a1**2
  a3squared = sim_a3**2
  a1squaredInv = 1.0/a1squared
  a3squaredInv = 1.0/a3squared

  ! Coefficients A1=A2, A3
  ! First, the limiting values
  if (abs(sim_e) < tiny(0.0)) then               ! sim_e=0, a sphere
     AA1ch = 2.0 / 3.0
     AA3ch = AA1ch
  else if (abs(sim_e - 1.0) < tiny(0.0) ) then   ! sim_e = 1, a disk?
     AA1ch = 0.0
     AA3ch = 2
  else                                            ! general case
     AA1ch = sqrt(1.0-sim_e**2)*asin(sim_e)/sim_e**3 - (1.0-sim_e**2)/sim_e**2
     AA3ch = 2.0/sim_e**2 - 2.0*sqrt(1.0-sim_e**2)*asin(sim_e)/sim_e**3
  end if
  Ich   = 2.0*AA1ch*a1squared + AA3ch*a3squared

  DD = sim_pi*sim_Newton*sim_density
  Mass = 4.0/3.0*sim_pi*a1squared*sim_a3*sim_density

  lo = tileDesc%limits(LOW,  :)
  hi = tileDesc%limits(HIGH, :)
  allocate(xLeft(lo(IAXIS):hi(IAXIS)))
  allocate(yLeft(lo(JAXIS):hi(JAXIS)))
  allocate(zLeft(lo(KAXIS):hi(KAXIS)))

  if (NDIM == 3) then  
     call Grid_getCellCoords(KAXIS, LEFT_EDGE, tileDesc%level, lo, hi, zLeft)
  endif
  if (NDIM >= 2) then    
     call Grid_getCellCoords(JAXIS, LEFT_EDGE, tileDesc%level, lo, hi, yLeft)
  endif
  call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level, lo, hi, xLeft)

  ! delta x is constant throughout each block
  call tileDesc%deltas(deltas)
  dx = deltas(IAXIS)
  dy = deltas(JAXIS)
  dz = deltas(KAXIS)

  !! Loop over all LEAF elements and generate the analytical solution
!  print *,'sim_nsubinv,nsubzones',sim_nsubinv,sim_nsubzones
  do       k = lo(KAXIS), hi(KAXIS)
     dzz = dz * sim_nsubinv
     do    j = lo(JAXIS), hi(JAXIS)
        dyy = dy * sim_nsubinv
        do i = lo(IAXIS), hi(IAXIS)
           dxx = dx * sim_nsubinv

           phiSum = 0.0

           ! volume average over subzones
           do kk = 0, sim_nsubzones-1
              zz    = zLeft(k) + (kk + 0.5)*dzz 
              zdist = (zz - sim_zctr) * K3D

              do jj = 0, sim_nsubzones-1
                 yy    = yLeft(j) + (jj + 0.5)*dyy
                 ydist = (yy - sim_yctr) * K2D

                 do ii = 0, sim_nsubzones-1
                    xx    = xLeft(i) + (ii + 0.5)*dxx
                    xdist = xx - sim_xctr

                    select case (sim_initGeometry)

                    case (POLAR)   ! 2d axisymmetric?

                       rFunc = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       r2 = xdist**2 + ydist**2

                    case (CYLINDRICAL)   ! 2d axisymmetric

                       rFunc = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       r2 = xdist**2
                       z2 = ydist**2

                    case (CARTESIAN)       ! 3d cartesian

                       r2 = xdist**2 + ydist**2
                       z2 = zdist**2
                       rFunc = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                    end select

                    if (rFunc <= 1.) then    ! inside the spheroid
                       phiSum = phiSum + DD*(Ich - AA1ch*r2 - AA3ch*z2) 
                    else                     ! outside the spheroid
                       lambda = 0.5d0*( (r2 + z2 - a1squared - a3squared) + &
                            &                        sqrt( (a1squared+a3squared-r2-z2)**2 -              &
                            &                        4.0*(a1squared*a3squared-r2*a3squared-z2*a1squared) ) )
                       Rinv = 1.0/sqrt(a3squared + lambda)
                       if (abs(sim_e) < tiny(0.0)) then   ! case for e=0  outside the sphere
                          ! H0 = 1, H1 = 1, H3 = 1
                          distTerm = 1.0 - (1.0/3.0)*Rinv**2*(r2 + z2)
                       else                               ! general case outside the spheroid
                          h = sim_a1 * sim_e * Rinv
                          hinv = 1.0/h                 
                          atan_h = atan(h)
                          H0 = atan_h * hinv
                          H1 = (atan_h - h/(1.0+h**2))*hinv**3
                          H3 = (h - atan_h)*2.0*hinv**3

                          distTerm = ( H0 - 0.5d0*Rinv**2*(r2*H1 + z2*H3) )
                       endif
                       phiSum = phiSum + 1.50*sim_Newton*Mass*Rinv*distTerm
                    endif

                 enddo
              enddo
           enddo

           call tileDesc%getDataPtr(solnData, CENTER)
           solnData(APOT_VAR, i, j, k) = -phiSum/(sim_nsubzones**3)   
           call tileDesc%releaseDataPtr(solnData, CENTER)
        enddo
     enddo
  enddo

  ! Cleanup
  deallocate(xLeft)
  deallocate(yLeft)
  deallocate(zLeft)

end subroutine sim_analytical
