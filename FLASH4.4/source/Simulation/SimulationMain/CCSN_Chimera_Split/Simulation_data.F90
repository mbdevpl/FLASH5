!!****if* source/Simulation/SimulationMain/CCSN/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: CCSN_Chimera
!!  
!! PARAMETERS
!!
!!
!! NOTES
!!  
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  character(len=256),save :: chimera_model_file
  character(len=256),save :: progenitor_model_file
  real, save :: sim_max_r, sim_r_inner, sim_rho_inner
  logical, save :: sim_do_quad
  integer, save :: sim_nquad

  integer, save :: nvar_stored
  integer, parameter :: n1d_max = 10000 ! Max number of lines a file can have
  integer, save :: n1d_total ! Actual number of lines, calculated after input
  real, save :: sim_smlrho, sim_smallt,sim_smallx
  real,save :: xzn(n1d_max), volxzn(n1d_max), dvolxzn(n1d_max)
  real,save :: model_1d(n1d_max,NUNK_VARS)
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save  :: sim_windVel, sim_massLoss, sim_velMult
  character (len=4), save :: unklabels(UNK_VARS_BEGIN:UNK_VARS_END)

  integer, save :: sim_meshMe

  logical,save :: useCool
  real, save :: coolTemp

  logical, save :: sim_restart, sim_burnUpdateEint
  real, save :: sim_pointMass, sim_expEner, sim_holeRadius

  real, save :: sim_Enu
  real, save :: sim_rhoOne, sim_rhoTwo, sim_rhoThree
  real, save :: sim_yOne, sim_yTwo, sim_yc, sim_yThree

  logical, save :: sim_usePnotT

  real, save :: sim_shockRadTot, sim_shockRadNum
  logical, save :: sim_postBounce
  real, save :: sim_bounceTime
  real, save :: sim_massAccRate, sim_massAccNum

  real, save :: sim_maxDens

end module Simulation_data
