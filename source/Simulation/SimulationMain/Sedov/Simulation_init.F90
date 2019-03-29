!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
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
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over initial explosion region)
!!  sim_minRhoInit     Density floor for initial condition
!!  sim_rInit          Radius of region into which explosion energy is dumped
!!                     initially, used only if tinitial <= 0.
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!  sim_forceCenterDerefine  Try to force low refinement level around explosion center?
!!  sim_centerRefineLevel    Desired refinement level at center (if "forcing")
!!  sim_derefineRadius       Radius of center region to force derefinement
!!  sim_profFileName   Name of file from which to read a 1D Sedov solution for the
!!                     initial condition
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use ut_generalInterface,         ONLY : ut_getFreeFileUnit

  implicit none

#include "Flash.h"
#include "constants.h"

  logical :: threadBlockListBuild, threadWithinBlockBuild  
  real    :: vctr

  call Driver_getMype(MESH_COMM,   sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, sim_globalNumProcs)

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
!!$  print*,'sim_rhoAmbient is',sim_rhoAmbient
  call RuntimeParameters_get('sim_expEnergy', sim_expEnergy)
  call RuntimeParameters_get('sim_rInit', sim_rInit)
  call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  call RuntimeParameters_get('sim_minRhoInit', sim_minRhoInit)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallT', sim_smallT)
  call RuntimeParameters_get('smallu', sim_smallu)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('tinitial',sim_tInitial)
  call RuntimeParameters_get('sim_forceCenterDerefine',sim_forceCenterDerefine)
  call RuntimeParameters_get('sim_centerRefineLevel',sim_centerRefineLevel)
  call RuntimeParameters_get('sim_derefineRadius',sim_derefineRadius)
  call RuntimeParameters_get('sim_profFileName',sim_profFileName)
  call RuntimeParameters_get('sim_bcSetBdryVar',sim_bcSetBdryVar)
  call RuntimeParameters_get('sim_earliestLSTime',sim_earliestLSTime)
  call RuntimeParameters_get('sim_latestLSTime',sim_latestLSTime)
  call RuntimeParameters_get('sim_smallestNormRadius',sim_smallestNormRadius)
  call RuntimeParameters_get('sim_largestNormRadius',sim_largestNormRadius)
  call RuntimeParameters_get('sim_oneLevelIntegralsOnly',sim_oneLevelIntegralsOnly)
  call RuntimeParameters_get('sim_integralsLevel',       sim_integralsLevel)

  sim_useProfileFromFile = .FALSE.
  if (sim_nProfile > 1) then
     sim_useProfileFromFile = (trim(sim_profFileName) .NE. "/dev/null")
  end if


  sim_inSubZones = 1./real(sim_nSubZones)
  sim_inSubzm1   = 1./real(max(sim_nSubZones-1,1))
  sim_inszd      = sim_inSubZones**NDIM
  
  !
  !  Calculate the initial volume and interior pressure.
  !
  if (NDIM .eq. 1) then
     vctr = 2. * sim_rInit
  elseif (NDIM .eq. 2) then
     vctr = sim_pi * sim_rInit**2
  else
     vctr = 4./3.*sim_pi*sim_rInit**3
  endif
  
  sim_pExp    = (sim_gamma-1.) * sim_expEnergy / vctr   ! only used if tinitial .LE. 0.


  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", sim_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", sim_threadWithinBlock)

  if (sim_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadBlockList = .false.
  end if
  if (sim_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadWithinBlock = .false.
  end if

  if (sim_useProfileFromFile) then
     call sim_readProfile()
     if (sim_tinitial > 0.0) call sim_scaleProfile(sim_tinitial)
  end if
  sim_analyticTime = -HUGE(1.0)
  sim_analyticGen  = -1
  if (sim_meshMe == MASTER_PE) then
     print*,'sim_rhoAmbient is',sim_rhoAmbient
  end if
  
  ! Open file for writing the numerical solution discretized onto the FLASH grid:
  sim_fileUnitOutNum = ut_getFreeFileUnit()
  open(unit=sim_fileUnitOutNum, status='SCRATCH')
!!$
!!$  ! Write the file header:
!!$  if(sim_globalME == MASTER_PE) then
!!$     write(sim_fileUnitOutNum,'(a10,7a15)') &
!!$          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint (spec.)'
!!$  end if
  
  ! Open file for writing the analytical solution discretized onto the FLASH grid:
  sim_fileUnitOutAna = ut_getFreeFileUnit()
!!$  open(unit=sim_fileUnitOutAna, file="sedSol-ana.out", form="formatted", position='append')
!!$
!!$  ! Write the file header:
!!$  if(sim_globalME == MASTER_PE) then
!!$     write(sim_fileUnitOutAna,'(a10,7a15)') &
!!$          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint density'
!!$  end if
  close(sim_fileUnitOutNum)

end subroutine Simulation_init
