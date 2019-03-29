!!****if* source/Simulation/SimulationMain/unitTest/Multipole/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the Maclaurin spheroid fluid data (density only) for a specified block
!!  and also sets the analytical gravitational solution by calling a separate routine.
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  blockID : the identification number of the block to be initialize
!!
!!***

subroutine Simulation_initBlock (SolnData, block)

  use Simulation_data, ONLY : sim_xctr,               &
                              sim_yctr,               &
                              sim_zctr,               &
                              sim_nsubinv,            &
                              sim_nsubzones,          &
                              sim_initGeometry,       &
                              sim_a3inv,              &
                              sim_a1inv,              &
                              sim_density,            &
                              sim_smallRho,           &
                              GRID_3DCARTESIAN,       &
                              GRID_3DCYLINDRICAL,     &
                              GRID_2DCYLINDRICAL,     &
                              GRID_1DSPHERICAL,       &
                              GRID_2DSPHERICAL

  use Grid_interface,  ONLY : Grid_getBlkIndexLimits, &
                              Grid_getCellCoords,     &
                              Grid_putRowData,        &
                              Grid_getDeltas
  use sim_interface, ONLY : sim_initBlockAnalytical
  use block_metadata, ONLY : block_metadata_t

  
  implicit none

#include "constants.h"
#include "Flash.h"

  type(block_metadata_t), intent(in) :: block
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: blockID, lev

  logical  :: gcell=.true.

  integer  :: sizeX, sizeY, sizeZ
  integer  :: i, j, k, ii, jj, kk
  integer  :: ibeg,iend,jbeg,jend,kbeg,kend

  real     :: dx, dy, dz
  real     :: xdist, ydist, zdist, Rdist, rFunc, r2, z2
  real     :: xx, yy, zz, dxx, dyy, dzz
  real     :: theta
  real     :: sum_rho

  integer, dimension(MDIM)      :: startingPos
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real,    dimension(MDIM)   :: deltas

  real, dimension(:), allocatable :: xLeft, yLeft, zLeft
  real, dimension(:), allocatable :: rho
!
!
!    ...Loop over all cells in current block and initialize the cells with
!       the needed density values. Work on a row-by-row basis. Break each cell
!       into nsubzones^ndim sub-zones and average the results to get the values
!       for the cell. This prevents a blocky spheroid shell.
!
!
!!$  call Grid_getBlkIndexLimits (blockID,    &
!!$                               blkLimits,  &
!!$                               blkLimitsGC )

  blockID=block%id
!!$  blkLimits=block%limits
!!$  blkLimitsGC=block%limitsGC
  blkLimits=block%locallimits
  blkLimitsGC=block%locallimitsGC
  lev=block%level
  select case (sim_initGeometry)
!
!
!   ...The 3D cartesian case.
!
!         Full consideration of subzones in all 3 dimensions necessary,
!         as square of radius (r2) and square of z coordinate (z2) for
!         evaluating the MacLaurin density vary with all 3 dimensions.
!
     
  case (GRID_3DCARTESIAN)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1
    sizeZ = blkLimitsGC (HIGH,KAXIS) - blkLimitsGC (LOW,KAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))
    allocate (zLeft (sizeZ))

    call Grid_getCellCoords (KAXIS, block, LEFT_EDGE, gcell, zLeft, sizeZ)
    call Grid_getCellCoords (JAXIS, block, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, block, LEFT_EDGE, gcell, xLeft, sizeX)

    call Grid_getDeltas     (lev, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)
    dz = deltas(KAXIS)

    kbeg = blkLimitsGC (LOW, KAXIS)
    kend = blkLimitsGC (HIGH,KAXIS)
    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    allocate (rho (sizeX))

    do k = kbeg,kend
       dzz = dz * sim_nsubinv
       do j = jbeg,jend
          dyy = dy * sim_nsubinv
          do i = ibeg,iend
             dxx = dx * sim_nsubinv

             sum_rho = 0.0

             do kk = 0, sim_nsubzones-1
                zz    = zLeft (k) + (kk + 0.5) * dzz 
                zdist = zz - sim_zctr
                do jj = 0, sim_nsubzones-1
                   yy    = yLeft (j) + (jj + 0.5) * dyy
                   ydist = yy - sim_yctr
                   do ii = 0, sim_nsubzones-1
                      xx    = xLeft (i) + (ii + 0.5) * dxx
                      xdist = xx - sim_xctr

                      r2    = xdist**2 + ydist**2
                      z2    = zdist**2

                      rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                      if (rFunc <= 1.) then                        ! inside the spheroid
                          sum_rho = sum_rho + sim_density
                      end if

                   enddo
                enddo
             enddo

             rho (i) = max (sim_smallRho + sim_smallRho , sum_rho * sim_nsubinv**3)

          enddo

          startingPos (1) = ibeg
          startingPos (2) = j
          startingPos (3) = k

          call Grid_putRowData (block, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
 
       enddo
    enddo

    deallocate (rho)
    deallocate (xLeft)
    deallocate (yLeft)
    deallocate (zLeft)
!
!
!   ...The 3D cylindrical case.
!
!         Subzones need to be considered only for the radial and z dimension.
!         The angular dimension does not need to be split into subzones (square
!         of radius r2 and square of z coordinate (z2) for evaluating the MacLaurin
!         potential vary only with radial and z dimension component of 3D cylindrical
!         coordinate system.
!
!
  case (GRID_3DCYLINDRICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, block, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, block, LEFT_EDGE, gcell, xLeft, sizeX)
    call Grid_getDeltas     (lev, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    kbeg = blkLimitsGC (LOW, KAXIS)
    kend = blkLimitsGC (HIGH,KAXIS)
    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    allocate (rho (sizeX))

    do j = jbeg,jend
       dyy = dy * sim_nsubinv           ! z component
       do i = ibeg,iend
          dxx = dx * sim_nsubinv        ! radial component

          sum_rho = 0.0

          do jj = 0, sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             zdist = yy - sim_yctr
             do ii = 0, sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr

                r2 = Rdist ** 2
                z2 = zdist ** 2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                        ! inside the spheroid
                    sum_rho = sum_rho + sim_density
                end if

             enddo
          enddo

          rho (i) = max (sim_smallRho + sim_smallRho , sum_rho * sim_nsubinv**2)

       enddo

       do k = kbeg,kend           ! angular component -> no subzones -> identical potentials

          startingPos (1) = ibeg
          startingPos (2) = j
          startingPos (3) = k

          call Grid_putRowData (block, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
 
       enddo

    enddo

    deallocate (rho)
    deallocate (xLeft)
    deallocate (yLeft)
!
!
!   ...The 2D cylindrical case.
!
!  
  case (GRID_2DCYLINDRICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, block, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, block, LEFT_EDGE, gcell, xLeft, sizeX)

    call Grid_getDeltas     (lev, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    allocate (rho (sizeX))

    do j = jbeg,jend
       dyy = dy * sim_nsubinv
       do i = ibeg,iend
          dxx = dx * sim_nsubinv

          sum_rho = 0.0

          do jj = 0, sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             zdist = yy - sim_yctr
             do ii = 0, sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr

                r2    = Rdist**2
                z2    = zdist**2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                        ! inside the spheroid
                    sum_rho = sum_rho + sim_density
                end if

             enddo
          enddo

          rho (i) = max (sim_smallRho + sim_smallRho , sum_rho * sim_nsubinv**2)

       enddo

       startingPos(1) = ibeg
       startingPos(2) = j

       call Grid_putRowData (block, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
 
    enddo

    deallocate (rho)
    deallocate (xLeft)
    deallocate (yLeft)
!
!
!   ...The 2D spherical case.
!
!  
  case (GRID_2DSPHERICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, block, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, block, LEFT_EDGE, gcell, xLeft, sizeX)

    call Grid_getDeltas     (lev, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    allocate (rho (sizeX))

    do j = jbeg,jend
       dyy = dy * sim_nsubinv
       do i = ibeg,iend
          dxx = dx * sim_nsubinv

          sum_rho = 0.0

          do jj = 0, sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             theta = yy - sim_yctr
             do ii = 0, sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr                     ! Rdist is x^2 + y^2 + z^2 here!

                zdist = Rdist * cos (theta)

                z2    = zdist**2                          ! this is z^2
                r2    = Rdist**2 - z2                     ! this is x^2 + y^2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                        ! inside the spheroid
                    sum_rho = sum_rho + sim_density
                end if

             enddo
          enddo

          rho (i) = max (sim_smallRho + sim_smallRho , sum_rho * sim_nsubinv**2)

       enddo

       startingPos(1) = ibeg
       startingPos(2) = j

       call Grid_putRowData (block, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
 
    enddo

    deallocate (rho)
    deallocate (xLeft)
    deallocate (yLeft)

  case default

    call Driver_abortFlash ('Multipole unit test: unsupported geometry (Simulation_initBlock)!')

  end select
!
!
!    ...Now calculate the analytical gravitational solution on this block.
!
!
  call sim_initBlockAnalytical (block)
!
!
!   ...Ready!
!
!
  return
end subroutine Simulation_initBlock
