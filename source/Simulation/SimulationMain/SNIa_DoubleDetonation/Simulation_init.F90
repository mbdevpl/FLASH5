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

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, only : Logfile_stampMessage
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry, Grid_getDomainBoundBox
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use Simulation_interface, ONLY : Simulation_mapStrToInt

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  character(len=256) :: initialWDFileName
  character(len=1024) :: header_line
  character(len=32), allocatable :: header_names(:)
  character(len=32) :: tmp_name
  character(len=4) :: unk_name
  integer, parameter :: iachar_tab = 9
  integer, parameter :: iachar_hash = iachar('#')
  integer :: istat, i, j, iunk, meshGeom, iachar_buffer, i1, i2, nvars
  real :: mass, dmass, vol, dvol, r1, r2, dr
  real, dimension(LOW:HIGH,MDIM) :: boundBox

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'initialWDFile', initialWDFileName)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  call RuntimeParameters_get('xhe4_fluff', sim_xhe4Fluff)
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
  call RuntimeParameters_get('xo16_fluff', sim_xo16Fluff)
  call RuntimeParameters_get('xni56_fluff', sim_xni56Fluff)

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

  call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile')
  open(unit=2,file=initialWDFileName,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile')

  ! parse the header to get isotope names
  read(2,"(a)") header_line

  ! replace '#' and tabs with spaces
  do i = 1, len_trim(header_line)
     iachar_buffer = iachar(header_line(i:i))
     if ( iachar_buffer == iachar_hash .or. iachar_buffer == iachar_tab ) header_line(i:i) = ' '
  end do

  ! remove leading/extra whitespace
  header_line = adjustl(header_line)
  i1 = 0
  i2 = len_trim(header_line)
  do i = 1, i2-1
     i1 = i1 + 1
     header_line(i1:i1) = header_line(i:i)
     if ( header_line(i:i+1) == '  ' ) i1 = i1 - 1
  end do
  i1 = i1 + 1
  header_line(i1:i1) = header_line(i2:i2)
  if ( i1 < i2 ) header_line(i1+1:) = ' '

  ! count the number of vars in header
  nvars = count( [(header_line(i:i),i=1,len_trim(header_line))] == ' ' ) + 1
  allocate(header_names(nvars))
  read(header_line,*) header_names

  ! determine number of species in profile and extract names from header, map to UNK vars
  sim_wd_nspec = nvars - 3
  allocate(sim_wd_spec_name(sim_wd_nspec))
  allocate(sim_wd_spec2unk(sim_wd_nspec))
  sim_wd_spec2unk(1:sim_wd_nspec) = NONEXISTENT
  sim_wd_unk2spec(SPECIES_BEGIN:SPECIES_END) = NONEXISTENT
  do i = 1, sim_wd_nspec
     tmp_name = adjustl(header_names(i+3))
     unk_name = tmp_name(1:4)

     call Simulation_mapStrToInt(unk_name,iunk,MAPBLOCK_UNK)
     if ( iunk /= NONEXISTENT ) then
        sim_wd_spec2unk(i) = iunk
        sim_wd_unk2spec(iunk) = i
     end if
  end do


  read(2,*) sim_wd_npnts

  allocate(sim_wd_rad_tab (0:sim_wd_npnts+1))
  allocate(sim_wd_dens_tab(1:sim_wd_npnts+1))
  allocate(sim_wd_temp_tab(1:sim_wd_npnts+1))
  allocate(sim_wd_spec_tab(1:sim_wd_npnts+1,sim_wd_nspec))
  allocate(sim_wd_vol_tab (0:sim_wd_npnts+1))
  allocate(sim_wd_mass_tab(0:sim_wd_npnts+1))
  mass = 0.0
  vol = 0.0
  sim_wd_rad_tab(0) = 0.0
  sim_wd_vol_tab(0) = 0.0
  sim_wd_mass_tab(0) = 0.0
  do i = 1, sim_wd_npnts
     read(2,*) sim_wd_rad_tab(i), sim_wd_dens_tab(i), sim_wd_temp_tab(i), (sim_wd_spec_tab(i,j),j=1,sim_wd_nspec)
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
  if ( sim_densFluff <= 0.0 ) sim_densFluff = sim_wd_dens_tab(sim_wd_npnts)
  if ( sim_tempFluff <= 0.0 ) sim_tempFluff = sim_wd_temp_tab(sim_wd_npnts)
  sim_wd_dens_tab(sim_wd_npnts+1) = sim_densFluff
  sim_wd_temp_tab(sim_wd_npnts+1) = sim_tempFluff
  if ( sim_xhe4Fluff > 0.0 .or. sim_xc12Fluff > 0.0 .or. sim_xo16Fluff > 0.0 .or. sim_xni56Fluff > 0.0 ) then
     sim_wd_spec_tab(sim_wd_npnts+1,:) = 0.0
     if ( sim_wd_unk2spec(HE4_SPEC) > 0 ) sim_wd_spec_tab(sim_wd_npnts+1,sim_wd_unk2spec(HE4_SPEC)) = sim_xhe4Fluff
     if ( sim_wd_unk2spec(C12_SPEC) > 0 ) sim_wd_spec_tab(sim_wd_npnts+1,sim_wd_unk2spec(C12_SPEC)) = sim_xc12Fluff
     if ( sim_wd_unk2spec(O16_SPEC) > 0 ) sim_wd_spec_tab(sim_wd_npnts+1,sim_wd_unk2spec(O16_SPEC)) = sim_xo16Fluff
     if ( sim_wd_unk2spec(NI56_SPEC) > 0 ) sim_wd_spec_tab(sim_wd_npnts+1,sim_wd_unk2spec(NI56_SPEC)) = sim_xni56Fluff
  else
     sim_wd_spec_tab(sim_wd_npnts+1,:) = sim_wd_spec_tab(sim_wd_npnts,:)
     if ( sim_wd_unk2spec(HE4_SPEC) > 0 ) sim_xhe4Fluff = sim_wd_spec_tab(sim_wd_npnts,sim_wd_unk2spec(HE4_SPEC))
     if ( sim_wd_unk2spec(C12_SPEC) > 0 ) sim_xc12Fluff = sim_wd_spec_tab(sim_wd_npnts,sim_wd_unk2spec(C12_SPEC))
     if ( sim_wd_unk2spec(O16_SPEC) > 0 ) sim_xo16Fluff = sim_wd_spec_tab(sim_wd_npnts,sim_wd_unk2spec(O16_SPEC))
     if ( sim_wd_unk2spec(NI56_SPEC) > 0 ) sim_xni56Fluff = sim_wd_spec_tab(sim_wd_npnts,sim_wd_unk2spec(NI56_SPEC))
  end if
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

  sim_wd_vol_tab(sim_wd_npnts+1)  = 4.0*PI/3.0 * sim_wd_rad_tab(sim_wd_npnts+1)**3
  sim_wd_mass_tab(sim_wd_npnts+1) = sim_wd_mass_tab(sim_wd_npnts) + &
      sim_wd_dens_tab(sim_wd_npnts+1) * ( sim_wd_vol_tab(sim_wd_npnts+1) - sim_wd_volume )

end subroutine Simulation_init
