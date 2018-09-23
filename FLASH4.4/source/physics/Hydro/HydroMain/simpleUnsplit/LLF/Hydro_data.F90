!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/HLL/Hydro_data
!!
!! NAME
!!
!!  Hydro_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Hydro_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data that are specific to the HLL unit implementation.
!!
!!***
 
 
Module Hydro_data

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  logical, save :: hy_fluxCorrect,        hy_charLimiting,     &
                   hy_shockDetectOn,      hy_hybridRiemannOnly,&
                   hy_useDiffuse,         hy_useViscosity,     &
                   hy_useConductivity,    hy_ContactSteepening,&
                   hy_EOSforRiemann,      hy_flattening,       &
                   hy_forceHydroLimit = .false.,               &
                   hy_killdivb = .false.,                      &
                   hy_energyFixSwitch = .false.,               &
                   hy_entropy,                                 &
                   hy_upwindTVD,          hy_use_avisc,        &
                   hy_gravConsv,          hy_useGravHalfUpdate,&
                   hy_useGravPotUpdate,   hy_useGravity,       &
                   hy_updateHydroFluxes,  hy_useHydro,         &
                   hy_use3dFullCTU,       hy_addThermalFlux,   &
                   hy_useHybridOrder,     hy_useVaryingCFL,    &
                   hy_conserveAngMom 

  integer, save :: hy_geometry, hy_limiter, hy_irenorm
  integer, save :: hy_unsplitEosMode = MODE_DENS_EI
  integer, save :: hy_EosModeAfter = MODE_DENS_EI
  integer, save :: hy_RiemannSolver,hy_RiemannSolverLoc, hy_entropyFixMethod
  integer, parameter :: hy_numXN = NSPECIES+NMASS_SCALARS

  character(len=MAX_STRING_LENGTH), save :: hy_limiter_str,      &
                                            hy_RiemannSolver_str,&
                                            hy_entropyFixMethod_str

  real, save :: hy_cfl
  real, save :: hy_cfl_original
  real, save :: hy_tiny=1.e-32
  real, save :: hy_Rconst
  real, save :: hy_eswitch
  real, save :: hy_smalldens
  real, save :: hy_smallpres
  real, save :: hy_smallu
  real, save :: hy_LimitedSlopeBeta
  real, save :: hy_maxDiff
  real, save :: hy_cvisc
  real, save :: hy_hybridOrderKappa

  ! System of units used
  character(4), save :: hy_units

  ! Everybody should know these!
  integer, save :: hy_meshNumProcs, hy_meshMe

  ! Constants for non-dimensionalization
  real, save :: hy_xref
  real, save :: hy_tref
  real, save :: hy_dref
  real, save :: hy_vref
  real, save :: hy_pref
  real, save :: hy_eref
  real, save :: hy_qref
  real, save :: hy_bref
  real, save :: hy_mref
  real, save :: hy_gref
  real, save :: hy_nref
  real, save :: hy_kref

  ! Order of Accuracy
  ! hy_order=1 ==> First order Godunov
  ! hy_order=2 ==> Muscl-Hancock
  ! hy_order=3 ==> PPM
  ! hy_order=5 ==> WENO: Note when using WENO, users should use 6 guardcells 
  !                      along with increased sizes of nxb, nyb, and nzb that are
  !                      larger than 2*NGUARD=12 at least.
  integer, save :: hy_order, hy_transOrder, hy_3Torder
  integer, save :: hy_gcMaskSize=NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(NUNK_VARS+NDIM*NFACE_VARS),save :: hy_gcMask

  integer, dimension(NFLUXES), save :: hy_fluxCorVars

!#ifdef FLASH_UHD_3T
  real,    save :: hy_avogadro, hy_qele
  real,    save :: hy_speedOfLight
!#endif

  logical, save :: hy_threadBlockList = .false.
  logical, save :: hy_threadTileList = .false.
  logical, save :: hy_threadWithinBlock = .false.


  integer, save :: hy_3TMode = HY3T_NONE

End Module Hydro_data
