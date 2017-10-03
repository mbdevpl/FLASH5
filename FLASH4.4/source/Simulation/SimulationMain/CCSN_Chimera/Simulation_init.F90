!!****if* source/Simulation/SimulationMain/CCSN/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initialization routine for CCSN problem.  Reads in 
!!  runtime parameters.  Reads 1D model file and maps
!!  variable names to FLASH unk names.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine Simulation_init()
  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getGeometry, Grid_getDomainBoundBox
  !use gr_mpoleData, ONLY : gr_point_mass=>point_mass

  use quadrature_module, only: quad_init
  use chimera_model_module
  use model_interp_module

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: i, j, k
  integer :: meshGeom

  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: mass_chim, vol_chim, mass_prog, vol_prog
  real :: r, theta, phi, x, y, z, rcyl, zcyl, gr_max_r, point_mass
  real :: rlo, rhi, dr, dvol, dvolr, domega, domega_exclude, tmp1, tmp2, tmp3
  integer :: irho_inner

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call Grid_getGeometry(meshGeom)
  call Grid_getDomainBoundBox(boundBox)

  call RuntimeParameters_get( 'chimera_model_file', chimera_model_file)
  call RuntimeParameters_get( 'progenitor_model_file', progenitor_model_file)

  call RuntimeParameters_get( 'smlrho',   sim_smlrho)
  call RuntimeParameters_get( 'smallt',   sim_smallt)
  call RuntimeParameters_get( 'smallx',   sim_smallx)
  call RuntimeParameters_get( 'restart', sim_restart)

  call RuntimeParameters_get( 'max_r', sim_max_r)
  call RuntimeParameters_get( 'r_inner', sim_r_inner)
  call RuntimeParameters_get( 'rho_inner', sim_rho_inner)
  call RuntimeParameters_get( 'do_quad', sim_do_quad)
  call RuntimeParameters_get( 'nquad', sim_nquad)

!  call RuntimeParameters_get('xhe4', sim_xhe4)
!  call RuntimeParameters_get('xc12', sim_xc12)
!  call RuntimeParameters_get('xo16', sim_xo16)


  ! load chimera model
  if ( len_trim(chimera_model_file) > 0 ) then
     call open_chimera_file(chimera_model_file)
     call read_chimera_file
  end if

  ! load progenitor model
  if ( len_trim(progenitor_model_file) > 0 ) then
    call read_progenitor_file(progenitor_model_file)
  end if

  ! set the maximum radius to use chimera data
  ! if r > sim_max_r, use progenitor data
  ! by default, use the inner-edge of outer-most chimera zone (excluding ghost zones)
  ! this obviates the need for any extrapolation of chimera data
  if ( sim_max_r <= 0.0 ) then
!    sim_max_r = (3.0 * volx_e_chim(imax_chim) )**(1.0/3.0)
     sim_max_r = x_e_chim(imax_chim)
  end if

  ! if sim_r_inner is explicitly defined at runtime, use that as the inner radius
  ! otherwise, use the radius corresponding to the inner edge of the first cell with
  ! density less than sim_rho_inner (default is 1.0e11 g/cm^3)
  if ( sim_r_inner <= 0.0 ) then
    irho_inner = locate_descending( sim_rho_inner, imax_chim, rhobar_c_chim )
    sim_r_inner = x_e_chim(irho_inner)
  end if

  ! set the inner radius to be the max of the user-defined value or the inner grid boundary
  sim_r_inner = max( boundBox(LOW,IAXIS), sim_r_inner )

  ! largest radius of FLASH grid
  if ( meshGeom == SPHERICAL ) then
     gr_max_r = boundBox(HIGH,IAXIS)
  else if ( meshGeom == CARTESIAN ) then
     gr_max_r = sqrt(sum(max(abs(boundBox(LOW,1:NDIM)),abs(boundBox(HIGH,1:NDIM)))**2))
  else if ( NDIM == 2 .and. meshGeom == CYLINDRICAL ) then
     gr_max_r = sqrt(max(abs(boundBox(LOW,JAXIS)),abs(boundBox(HIGH,JAXIS)))**2+boundBox(HIGH,IAXIS)**2)
  end if

  mass_chim = 0.0
  vol_chim = 0.0
  if ( NDIM == 1 .and. meshGeom == SPHERICAL ) then
     ! integrate mass from chimera data on FLASH grid (1D spherical only)
     do i = imin_chim, imax_chim
        if ( x_e_chim(i) < min(boundBox(HIGH,IAXIS),sim_max_r) ) then
           if ( x_e_chim(i+1) <= min(boundBox(HIGH,IAXIS),sim_max_r) ) then
              mass_chim = mass_chim + sum(dmass_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim))
              vol_chim = vol_chim + sum(dvol_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim))
           else
              dr = x_e_chim(i+1) - min(boundBox(HIGH,IAXIS),sim_max_r)
              dvolr = dr * ( x_e_chim(i) * min(boundBox(HIGH,IAXIS),sim_max_r) + dr * dr * (1.0/3.0) )
              dvol  = 4.0 * PI * ( dvolx_e_chim(i) - dvolr )

              vol_chim = vol_chim + dvol
              mass_chim = mass_chim + ( dvolx_e_chim(i) - dvolr ) * &
              &                       sum( rho_c_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) * &
              &                            domega_chim(jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
           end if
        end if
     end do
  else if ( NDIM == 2 .and. meshGeom == CYLINDRICAL ) then
     ! integrate mass from chimera data on FLASH grid (2D cylindrical only)
     do j = jmin_chim, jmax_chim
        do i = imin_chim, imax_chim
           if ( x_e_chim(i) < sim_max_r ) then
              rcyl = x_c_chim(i) * sin( y_c_chim(j) )
              zcyl = x_c_chim(i) * cos( y_c_chim(j) )
              if ( rcyl > boundBox(LOW,IAXIS) .and. rcyl < boundBox(HIGH,IAXIS) .and. &
                   zcyl > boundBox(LOW,JAXIS) .and. zcyl < boundBox(HIGH,JAXIS) ) then
                 if ( x_e_chim(i+1) <= sim_max_r ) then
                    mass_chim = mass_chim + sum(dmass_e_chim(i,j,kmin_chim:kmax_chim))
                    vol_chim = vol_chim + sum(dvol_e_chim(i,j,kmin_chim:kmax_chim))
                 else
                    dr = x_e_chim(i+1) - sim_max_r
                    dvolr = dr * ( x_e_chim(i) * sim_max_r + dr * dr * (1.0/3.0) )
                    dvol  = sum(domega_chim(j,kmin_chim:kmax_chim)) * ( dvolx_e_chim(i) - dvolr )

                    vol_chim = vol_chim + dvol
                    mass_chim = mass_chim + ( dvolx_e_chim(i) - dvolr ) * &
                    &                       sum( rho_c_chim(i,j,kmin_chim:kmax_chim) * &
                    &                            domega_chim(j,kmin_chim:kmax_chim) )
                 end if
              end if
           end if
        end do
     end do
  else if ( NDIM == 3 .and. meshGeom == CARTESIAN ) then
     ! integrate mass from chimera data on FLASH grid (3D cartesian only)
     do k = kmin_chim, kmax_chim
        phi = z_c_chim(k)
        do j = jmin_chim, jmax_chim
           theta = y_c_chim(j)
           do i = imin_chim, imax_chim
              r = x_c_chim(i)
              rlo = x_e_chim(i)
              rhi = x_e_chim(i+1)
              ! only use chimera data if zone is inside sim_max_r
              if ( rlo < sim_max_r ) then
                 x = r * sin( theta ) * cos( phi )
                 y = r * sin( theta ) * sin( phi )
                 z = r * cos( theta )
                 ! only include in sum if the zone is on FLASH grid
                 if ( x > boundBox(LOW,IAXIS) .and. x < boundBox(HIGH,IAXIS) .and. &
                 &    y > boundBox(LOW,JAXIS) .and. y < boundBox(HIGH,JAXIS) .and. &
                 &    z > boundBox(LOW,KAXIS) .and. z < boundBox(HIGH,KAXIS) ) then
                    if ( rhi > sim_max_r ) then
                       dr = rhi - sim_max_r
                       dvolr = dr * ( rlo * sim_max_r + dr * dr * (1.0/3.0) )
                       dvol  = domega_chim(j,k) * ( dvolx_e_chim(i) - dvolr )

                       vol_chim = vol_chim + dvol
                       mass_chim = mass_chim + rho_c_chim(i,j,k) * dvol
                    else
                       mass_chim = mass_chim + dmass_e_chim(i,j,k)
                       vol_chim = vol_chim + dvol_e_chim(i,j,k)
                    end if
                 end if
              end if
           end do
        end do
     end do
  end if

  mass_prog = 0.0
  vol_prog = 0.0
  ! integrate mass from progenitor data on FLASH grid
  do i = 1, n1d_total
     rlo = xzn(i)
     rhi = xzn(i+1)
     ! this includes progenitor data on FLASH grid but not covered by chimera data
     if ( rhi > sim_max_r .and. rlo < gr_max_r ) then
        r = 0.5 * ( rlo + rhi )
        dvolr = 0.0
        if ( rhi > gr_max_r ) then
           dr = rhi - gr_max_r
           dvolr = dvolr + dr * ( gr_max_r * rhi + dr * dr * (1.0/3.0) )
        end if
        if ( rlo < sim_max_r ) then
           dr = sim_max_r - rlo
           dvolr = dvolr + dr * ( rlo * sim_max_r + dr * dr * (1.0/3.0) )
        end if

        if ( NDIM == 1 .and. meshGeom == SPHERICAL ) then
           domega = 4.0 * PI
        else if ( NDIM == 2 .and. meshGeom == CYLINDRICAL ) then
           tmp1 = min( 1.0, max( -1.0, boundBox(HIGH,JAXIS)/r ) )
           tmp2 = min( 1.0, max( -1.0, boundBox(LOW,JAXIS)/r ) )
           tmp3 = 2.0 * sqrt( 1.0 - min( 1.0, boundBox(HIGH,IAXIS)/r )**2 )
           domega = 2.0 * PI * ( tmp1 - tmp2 - tmp3 )
           domega = min( 4.0 * PI, max( 0.0, domega ) )
        else if ( NDIM == 3 .and. meshGeom == CARTESIAN ) then
           ! TODO: calculate the solid angle of shell exterior to FLASH grid
           domega_exclude = 0.0
           domega = 4.0 * PI - domega_exclude
        end if
        dvol = domega * ( dvolxzn(i) - dvolr )

        vol_prog = vol_prog + dvol
        mass_prog = mass_prog + dvol * model_1d(i,DENS_VAR)
     end if
  end do

  if (sim_meshMe == MASTER_PE) then
     write(*,'(a,2es23.15)') 'total mass   (chimera, progenitor)  =', mass_chim, mass_prog
     write(*,'(a,es23.15)')  'total mass   (chimera + progenitor) =', mass_chim + mass_prog

     write(*,'(a,2es23.15)') 'total volume (chimera, progenitor)  =', vol_chim, vol_prog
     write(*,'(a,es23.15)')  'total volume (chimera + progenitor) =', vol_chim + vol_prog
     if ( NDIM == 1 .and. meshGeom == SPHERICAL ) then
        dr = boundBox(HIGH,IAXIS) - boundBox(LOW,IAXIS)
        write(*,'(a,es23.15)') 'total volume (FLASH)         =', &
        &                       4.0 * PI * dr * ( boundBox(LOW,IAXIS) * boundBox(HIGH,IAXIS) + dr * dr * (1.0/3.0) )
     else if ( NDIM == 2 .and. meshGeom == CYLINDRICAL ) then
        write(*,'(a,es23.15)') 'total volume (FLASH)         =', &
        &                       PI * ( boundBox(HIGH,IAXIS)**2 - boundBox(LOW,IAXIS)**2 ) &
        &                          * ( boundBox(HIGH,JAXIS)    - boundBox(LOW,JAXIS) )
     else if ( NDIM == 3 .and. meshGeom == CARTESIAN ) then
        write(*,'(a,es23.15)') 'total volume (FLASH)         =', product( boundBox(HIGH,:) - boundBox(LOW,:) )
     end if
  end if

  ! set up quadrature weights and abscissae
  if ( sim_do_quad ) call quad_init( sim_nquad )

  point_mass = 0.0
  do i = 1, imax_chim+1
     if ( x_e_chim(i) < sim_r_inner ) then
        if ( x_e_chim(i+1) <= sim_r_inner ) then
           point_mass = point_mass + sum( dmass_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
        else
           dr = x_e_chim(i+1) - sim_r_inner
           dvolr = dr * ( sim_r_inner * x_e_chim(i+1) + dr * dr * (1.0/3.0) )
           point_mass = point_mass + ( dvolx_e_chim(i) - dvolr ) &
           &                         * sum( rho_c_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) &
           &                                * domega_chim(jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
        end if
     end if
  end do

  if ( sim_meshMe == MASTER_PE ) then
     write(*,*) 'point_mass=',point_mass
     write(*,*) 'adding this value to user-supplied value'
  end if
  !gr_point_mass = gr_point_mass + point_mass
  sim_pointMass = point_mass
  contains

!!! THIS IS OLD READ IN OF KEPLER FILE
!!!----------------------------------------------------------------------------
     subroutine read_progenitor_file(progenitor_model_file)

        use Simulation_interface, ONLY : Simulation_mapIntToStr

        implicit none

        ! input variables
        character (*), intent(in) :: progenitor_model_file

        ! local variables
        integer, parameter :: max_stored_vars = 30
        real :: var_temp(max_stored_vars)
        integer :: ipos, NUNK_VARS_stored
        integer :: var_key(NUNK_VARS)
        character (len=4) :: var_labels(max_stored_vars)
        character (len=256) :: current_line
        integer :: prog_unit
        real :: dr
        integer :: i, j, ierr

        do i = UNK_VARS_BEGIN, UNK_VARS_END
           call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
           call makeLowercase(unklabels(i))
        enddo

        ! open the file and read in the header 
        open (newunit=prog_unit,file=trim(progenitor_model_file),status='old')
        read (prog_unit,'(a256)') current_line

        if (sim_meshMe == MASTER_PE) then
           print *, 'file opened'
        end if

        ! read in the number of variables line
        read (prog_unit,'(a256)') current_line
        ipos = index(current_line,'=') + 1
        read (current_line(ipos:),*) nvar_stored
        if (sim_meshMe == MASTER_PE) then
           print *,"read nvar_stored", nvar_stored
        end if

        if (NUNK_VARS /= nvar_stored .AND. sim_meshMe == MASTER_PE) then
           print *, ' '
           print *, 'Warning: the number of variables stored in the'
           print *, 'input file is different than the number of'
           print *, 'variables in the current version of FLASH.'
           print *, ' '
           print *, 'The variables in the file that are also defined'
           print *, 'in FLASH will be read in.  Any missing variables'
           print *, 'will be initialized to zero'
           print *, ' '
        endif

        if (sim_meshMe == MASTER_PE) then
           print *, "Vaiables in file:"
        endif

        do i = 1, nvar_stored
           read (prog_unit,'(a4)') var_labels(i)
           if (sim_meshMe == MASTER_PE) then
              print *, var_labels(i)
           end if
           call makeLowercase(var_labels(i))
        enddo

        do j = 1, NUNK_VARS
           var_key(j) = NONEXISTENT
           do i = 1, nvar_stored
              if (unklabels(j) == var_labels(i)) then
                 var_key(j) = i
              endif
           enddo
           if (var_key(j) == NONEXISTENT) then
              if(sim_meshMe == MASTER_PE) then
                 print *, 'Warning, variable: ', unklabels(j), ' not found in the input file.'
                 print *, 'initializing ', unklabels(j), ' to 0'
                 print *, ' '
              endif
           endif
        enddo

        do i = 1, n1d_max
           read(prog_unit,*,iostat=ierr) xzn(i), (var_temp(j),j=1,nvar_stored)
           if ( ierr /= 0 ) exit
           ! put these in order, so model1d_var always contains the same variables in the same spots
           do j = 1, NUNK_VARS
              if (var_key(j) /= NONEXISTENT) then
                 model_1d(i,j) = var_temp(var_key(j))
              else
                 model_1d(i,j) = 0.0
              endif
           enddo
        enddo
        close(prog_unit)

        n1d_total = 0
        do while (xzn(n1d_total+1) /= 0.0)
           n1d_total = n1d_total + 1
        enddo

        if (sim_meshMe == MASTER_PE) then
           print *, 'file read completed'
           print *, n1d_total, 'points read in'
        endif

        volxzn(1) = (1.0/3.0) * xzn(1)**3
        do j = 1, n1d_total
           dr          = xzn(j+1) - xzn(j)
           dvolxzn(j)  = dr * ( xzn(j) * xzn(j+1) + dr * dr * (1.0/3.0) )
           volxzn(j+1) = volxzn(j) + dvolxzn(j)
        end do

  end subroutine read_progenitor_file
!!!----------------------------------------------------------------------------

end subroutine Simulation_init
