! Dean M. Townley 2009
!
! This is the static module data for the Simulation Unit for the SNIa_ddt setup
! Note that some of these are allocatable, and are allocated by Simulation_init()

module Simulation_data
#include "Flash.h"
#include "Eos.h"

  real,allocatable,dimension(:),save :: sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_c12_tab, sim_wd_ne22_tab
  real,allocatable,dimension(:),save :: sim_wd_rad_tab, sim_wd_vol_tab, sim_wd_mass_tab
  real, save :: sim_wd_radius, sim_wd_volume, sim_wd_mass, sim_wd_dr, sim_wd_dr_inv
  integer, save :: sim_wd_npnts

  logical, save :: sim_useShell
  real, save :: sim_radShellMin, sim_radShellMax, sim_thtShellMin, sim_thtShellMax
  real, save :: sim_xhe4Shell, sim_xc12Shell, sim_xni56Shell
  real, save :: sim_densShellMult, sim_tempShellMult
  real, save :: sim_densShell, sim_tempShell

  real, save :: sim_densFluff, sim_tempFluff, sim_xc12Fluff, sim_xne22Fluff

  logical, save :: sim_ignite
  real, save :: sim_ignX, sim_ignY, sim_ignZ
  real, save :: sim_ignRInner, sim_ignROuter
  real, save :: sim_ignTInner, sim_ignTOuter

  real, save :: sim_smallrho, sim_smallt, sim_smallp, sim_smalle
  real, save :: sim_smallx

  integer, save :: sim_globalMe

end module Simulation_data
