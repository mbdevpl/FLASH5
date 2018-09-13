!!
!! Dean M. Townsley 2009
!!
!! Initialization of Simulation Unit for SNIa_ddt.
!! See source/Simulation/Simulation_init.F90 for API spec and notes.
!! 
!! This routine does the following:
!! Read in initial WD and store in Simulation Unit static module data area.
!! Retrieve refinement parameters from parameter file datastructures.
!! Retrieve parameters for initial condition and do further setup if necessary.
!!
!! If a randomized initial condition is chosen, the coefficients,
!! which are global, are calculated here.  They are simply calculated on
!! every processor because they are fairly cheap.
!! This initial condition is the randomized multipole used to generate
!! the sample in Townsley etal (2009, ApJ, in press, ArXiv:0906.4384).

subroutine Simulation_init(myPE)

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, only : Logfile_stampMessage
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry, Grid_getDomainBoundBox
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  integer,intent(IN) :: myPE

  character(len=256) :: initialWDFileName
  integer :: istat, i, meshGeom
  real :: mass, dmass, vol, dvol, r1, r2, dr
  real, dimension(LOW:HIGH,MDIM) :: boundBox

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'initialWDFile', initialWDFileName)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
  call RuntimeParameters_get('xne22_fluff', sim_xne22Fluff)

  call RuntimeParameters_get('ignite', sim_ignite)
  call RuntimeParameters_get('x_match', sim_ignX)
  call RuntimeParameters_get('y_match', sim_ignY)
  call RuntimeParameters_get('z_match', sim_ignZ)
  call RuntimeParameters_get('r_match_inner', sim_ignRInner)
  call RuntimeParameters_get('r_match_outer', sim_ignROuter)
  call RuntimeParameters_get('t_ignite_inner', sim_ignTInner)
  call RuntimeParameters_get('t_ignite_outer', sim_ignTOuter)

  call RuntimeParameters_get('useShell', sim_useShell)
  call RuntimeParameters_get('radShellMin', sim_radShellMin)
  call RuntimeParameters_get('radShellMax', sim_radShellMax)
  call RuntimeParameters_get('thtShellMin', sim_thtShellMin)
  sim_thtShellMin = sim_thtShellMin * PI / 180.0
  call RuntimeParameters_get('thtShellMax', sim_thtShellMax)
  sim_thtShellMax = sim_thtShellMax * PI / 180.0
  call RuntimeParameters_get('densShellMult', sim_densShellMult)
  call RuntimeParameters_get('tempShellMult', sim_tempShellMult)
  call RuntimeParameters_get('densShell', sim_densShell)
  call RuntimeParameters_get('tempShell', sim_tempShell)

  call RuntimeParameters_get('xhe4_shell', sim_xhe4Shell)
  call RuntimeParameters_get('xc12_shell', sim_xc12Shell)
  call RuntimeParameters_get('xni56_shell', sim_xni56Shell)
  
  call RuntimeParameters_get('smlrho', sim_smallrho)
  call RuntimeParameters_get('smallt', sim_smallt)
  call RuntimeParameters_get('smallp', sim_smallp)
  call RuntimeParameters_get('smalle', sim_smalle)
  call RuntimeParameters_get('smallx', sim_smallx)

  call Driver_getMype(GLOBAL_COMM,sim_globalMe)

  !--------------------------------------------------------
  !  read in 1d initial wd profile
  !--------------------------------------------------------

!  call Logfile_stampMessage(myPE, '[Simulation_init] Reading initial 1-d WD profile')
  call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile')
  open(unit=2,file=initialWDFileName,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile')

  ! eat header
  read(2,*)
  read(2,*) sim_wd_npnts

  allocate(sim_wd_rad_tab (0:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_dens_tab(1:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_temp_tab(1:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_c12_tab (1:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_ne22_tab(1:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_vol_tab (0:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_mass_tab(0:sim_wd_npnts+1),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  mass = 0.0
  vol = 0.0
  sim_wd_rad_tab(0) = 0.0
  sim_wd_vol_tab(0) = 0.0
  sim_wd_mass_tab(0) = 0.0
  do i = 1, sim_wd_npnts
     read(2,*) sim_wd_rad_tab(i), sim_wd_dens_tab(i), sim_wd_temp_tab(i), sim_wd_c12_tab(i), sim_wd_ne22_tab(i)
     r1 = sim_wd_rad_tab(i-1)
     r2 = sim_wd_rad_tab(i)
     dr = r2 - r1
     dvol = dr * ( 3.0*r1*r2 + dr*dr ) * 4.0*PI/3.0
     dmass = sim_wd_dens_tab(i) * dvol
     vol = vol + dvol
     mass = mass + dmass
     sim_wd_vol_tab(i) = vol
     sim_wd_mass_tab(i) = mass
  enddo
  close(2)
  sim_wd_radius = sim_wd_rad_tab(sim_wd_npnts)
  sim_wd_volume = sim_wd_vol_tab(sim_wd_npnts)
  sim_wd_mass = sim_wd_mass_tab(sim_wd_npnts)
  sim_wd_dr = dr
  sim_wd_dr_inv = 1.0 / dr

  call Grid_getGeometry(meshGeom)
  call Grid_getDomainBoundBox(boundBox)
  if ( meshGeom == SPHERICAL ) then
     sim_wd_rad_tab(sim_wd_npnts+1)  = boundBox(HIGH,IAXIS)
  else if ( meshGeom == CYLINDRICAL ) then
     sim_wd_rad_tab(sim_wd_npnts+1)  = sqrt( sum( maxval(boundBox(LOW:HIGH,IAXIS:JAXIS)**2,dim=1) ) )
  else if ( meshGeom == CARTESIAN ) then
     sim_wd_rad_tab(sim_wd_npnts+1)  = sqrt( sum( maxval(boundBox**2,dim=1) ) )
  else
     call Driver_abortFlash("Geometry not supported")
  end if
  sim_wd_dens_tab(sim_wd_npnts+1) = sim_densFluff
  sim_wd_temp_tab(sim_wd_npnts+1) = sim_tempFluff
  sim_wd_c12_tab(sim_wd_npnts+1)  = sim_xc12Fluff
  sim_wd_ne22_tab(sim_wd_npnts+1) = sim_xne22Fluff
  sim_wd_vol_tab(sim_wd_npnts+1)  = 4.0*PI/3.0 * sim_wd_rad_tab(sim_wd_npnts+1)**3
  sim_wd_mass_tab(sim_wd_npnts+1) = sim_wd_mass_tab(sim_wd_npnts) + &
                                                            sim_densFluff * ( sim_wd_vol_tab(sim_wd_npnts+1) - sim_wd_volume )

end subroutine Simulation_init
