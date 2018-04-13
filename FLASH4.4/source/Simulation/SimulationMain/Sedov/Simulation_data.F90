!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
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
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over initial explosion region)
!!  sim_rInit          Radius of region into which explosion energy is dumped
!!                     initially, used only if tinitial <= 0.
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  real,parameter :: sim_pi = PI

  !! *** Runtime Parameters *** !!

  real, save    :: sim_pAmbient, sim_rhoAmbient, sim_expEnergy, sim_rInit
  real, save    :: sim_gamma, sim_xCenter, sim_yCenter, sim_zCenter
  real, save    :: sim_smallX, sim_smallRho, sim_smallP, sim_minRhoInit
  real, save    :: sim_smallT
  real, save    :: sim_smallu
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save :: sim_nSubZones
  real, save    :: sim_tInitial
  logical, save :: sim_forceCenterDerefine
  integer, save :: sim_centerRefineLevel
  real, save    :: sim_derefineRadius
  logical, save :: sim_bcSetBdryVar
  logical, save ::  sim_oneLevelIntegralsOnly
  integer,save  ::  sim_integralsLevel = 0

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter                  :: sim_nProfile = 10000
  real   , save                       :: sim_inSubZones, sim_inSubzm1
  real   , save                       :: sim_inszd
  real, dimension(4,sim_nProfile+1), save :: sim_profileInitial
  real, dimension(  sim_nProfile+1), save :: sim_rProf, sim_rhoProf, sim_pProf
  real, dimension(  sim_nProfile+1), save :: sim_vProf
  real, save                              :: sim_pExp
  logical, save                           :: sim_useProfileFromFile = .FALSE.
  character(len=MAX_STRING_LENGTH),  save ::          sim_profFileName
  logical,save                        ::  sim_profileIsScaled = .FALSE.
  real ,save                          ::  sim_profileScaledTime=0.0
  real ,save                          ::  sim_analyticTime
  integer,save                        ::  sim_analyticGen

  integer, save :: sim_meshMe, sim_globalMe
  integer, save ::             sim_globalNumProcs
  logical, save :: sim_threadBlockList = .false.
  logical, save :: sim_threadWithinBlock = .false.

  integer, save                       :: sim_fileUnitOutNum = 91, sim_fileUnitOutAna = 92 !initial values not really used
  character(len=80)                   :: sim_ffNum, sim_ffAna

#ifdef MAGP_VAR
  ! Cf. IO_writeIntegralQuantities, should be nGlobalSumProp+2*nRelErrSumProp+nGlobaExtProp
  integer,parameter :: nLargestMaxSummary = 20 + 5 + 3 + 15
#else
  integer,parameter :: nLargestMaxSummary = 19 + 5 + 3 + 15
#endif
  logical, save :: sim_testInitialized
  real, save    :: sim_testFirstVals(0:nLargestMaxSummary)
  real, save    :: sim_testLastVals(0:nLargestMaxSummary)
  real, save    :: sim_testLargestVals(0:nLargestMaxSummary)
  real, save    :: sim_testLargestWhen(0:nLargestMaxSummary)

  real, save    :: sim_earliestLSTime, sim_latestLSTime
  real, save    :: sim_smallestNormRadius, sim_largestNormRadius
end module Simulation_data
