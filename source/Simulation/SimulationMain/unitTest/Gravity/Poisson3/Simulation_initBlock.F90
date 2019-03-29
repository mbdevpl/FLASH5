!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: tileDesc) 
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified tileDesc.  This version sets up the Maclaurin spheroid problem.
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  tileDesc    --   The number of the tile to initialize
!!  
!!
!!***

!!REORDER(4): solnData
subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data, ONLY: sim_xctr, sim_yctr, sim_zctr, &
       &   sim_nsubinv, sim_nsubzones, sim_initGeometry, &
       &   sim_a3inv, sim_a1inv, sim_Pconst, sim_Omega2, sim_density, &
       &   sim_smallRho, sim_smallP, sim_gamma, sim_smallE
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "constants.h"
#include "Flash.h"

  real,pointer, dimension(:,:,:,:) :: solnData
  type(Grid_tile_t), intent(in)  :: tileDesc

  integer, dimension(MDIM) :: lo, hi
  real, dimension(MDIM)      :: deltas
  real     :: dx, dy, dz
  integer  :: i, j, k, ii, jj, kk
  real     :: xdist, ydist, zdist, dist2, rxy, rxyz2, rinv
  real     :: xx, yy, zz, dxx, dyy, dzz, vxfac, vyfac, vzfac
  real     :: sum_rho, sum_p, sum_vx, sum_vy, sum_vz, pres, vel

  real, dimension(:), allocatable :: xLeft, yLeft, zLeft
  real :: vx, vy, vz, p, rho, ek, ei


  ! Get the coordinate information for the current tile

  lo=tileDesc%limits(LOW,:)
  hi=tileDesc%limits(HIGH,:)
  
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

  ! Set initial conditions in each zone
  do       k = lo(KAXIS), hi(KAXIS)
     dzz = dz * sim_nsubinv
     do    j = lo(JAXIS), hi(JAXIS)
        dyy = dy * sim_nsubinv
        do i = lo(IAXIS), hi(IAXIS)
           dxx = dx * sim_nsubinv

           sum_rho = 0.0
           sum_p   = 0.0
           sum_vx  = 0.0
           sum_vy  = 0.0
           sum_vz  = 0.0

           ! Break the zone into nsubzones^ndim sub-zones and average the
           ! results to get the values for the zone.  This prevents a blocky
           ! spheroid shell

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

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       rxy   = xdist
                       rxyz2 = dist2
                       rinv  = 1./sqrt(xdist**2 + ydist**2)
                       vxfac = 0.0
                       vyfac = 0.0
                       vzfac = 0.0

                    case (CYLINDRICAL)   ! 2d axisymmetric

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       rxy   = xdist
                       rxyz2 = dist2
                       rinv  = 1./sqrt(xdist**2 + ydist**2)
                       vxfac = 0.0
                       vyfac = 0.0
                       vzfac = 0.0

                    case (CARTESIAN)       ! 3d cartesian

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rxy   = sqrt(xdist**2 + ydist**2)
                       rxyz2 = (rxy*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rinv  = 1./sqrt(xdist**2 + ydist**2 + zdist**2)
                       vxfac = -ydist*rinv
                       vyfac = xdist*rinv
                       vzfac = 0.0

                    end select

                    if (dist2 <= 1.) then    ! inside the spheroid
                       pres    = sim_Pconst * (1.0 - rxyz2)
                       vel     = rxy * sim_Omega2
                       sum_rho = sum_rho + sim_density
                       sum_p   = sum_p   + pres
                       sum_vx  = sum_vx  + vel*vxfac
                       sum_vy  = sum_vy  + vel*vyfac
                       sum_vz  = sum_vz  + vel*vzfac
                    else                     ! outside the spheroid
                       sum_rho = sum_rho + sim_smallRho*10
                       sum_p   = sum_p   + sim_smallP*10
                    endif

                 enddo
              enddo
           enddo

           rho = max(sum_rho * sim_nsubinv**3, sim_smallRho)
           p   = max(sum_p * sim_nsubinv**3, sim_smallP)
           vx  = sum_vx*sim_nsubinv**3
           vy  = sum_vy*sim_nsubinv**3
           vz  = sum_vz*sim_nsubinv**3
           ei  = max(p/(rho*(sim_gamma-1.0)), sim_smallE)
           ek  = 0.5*(vx**2+vy**2+vz**2)
           
           solnData(DENS_VAR, i, j, k) = rho
           solnData(PRES_VAR, i, j, k) = p
           solnData(VELX_VAR, i, j, k) = vx
           solnData(VELY_VAR, i, j, k) = vy
           solnData(VELZ_VAR, i, j, k) = vz
           solnData(GAME_VAR, i, j, k) = sim_gamma
           solnData(GAMC_VAR, i, j, k) = sim_gamma 
           solnData(EINT_VAR, i, j, k) = ei 
           solnData(ENER_VAR, i, j, k) = ei + ek

        enddo
     enddo
  enddo
  deallocate(xLeft)
  deallocate(yLeft)
  deallocate(zLeft)

  ! Now calculate the analytical solution on this tile
  
  call sim_analytical(tileDesc)

  return
end subroutine Simulation_initBlock

