!!****if* source/Grid/GridMain/AMR/Amrex/Grid_data
!!
!! NAME
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use Grid_data
!!
!! DESCRIPTION 
!!  
!!  Defining variables for storing AMReX-related data.
!!   
!!***

#include "Flash.h"
#include "constants.h"

Module Grid_data
  use iso_c_binding, ONLY : c_ptr

  implicit none

  integer, save :: gr_dirGeom(MDIM)
  real,    save :: gr_smalle, gr_smallrho

  logical, save :: gr_convertToConsvdForMeshCalls
  logical ,save :: gr_convertToConsvdInMeshInterp

  !!!!! NEEDED BY Grid_limitAbundance.F90
  real,    save :: gr_smallx

  !!!!! NEEDED BY gr_bcApplyToAllBlks.F90
  integer, save :: gr_numDataStruct
  integer, save :: gr_gridDataStruct(NDATATYPES)
  integer, save :: gr_gridDataStructSize(NDATATYPES)
  
  !!!!!  NEEDED BY gr_setDataStructInfo.F90
  logical, save :: gr_bcEnableApplyMixedGds

  ! Maintain local copies of AMReX-controlled data for optimization
  integer, save :: gr_lRefineMax
  integer, save :: gr_maxRefine
  integer, save :: gr_geometry
  integer, save :: gr_nrefs
  real,    save :: gr_minCellSize
  real,    save :: gr_minCellSizes(MDIM)
  logical, save :: gr_allPeriodic

  ! Local copies that stores BC information for AMReX callbacks.
  ! These variables should only be used by the AMReX callbacks.
  integer, target, save :: lo_bc_amrex(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)
  integer, target, save :: hi_bc_amrex(NDIM, UNK_VARS_BEGIN:UNK_VARS_END)

  ! These are historical.
  ! Within the AMReX implementation, the number of guardcells
  ! is set to NGUARD for all directions.  Code in the Amrex 
  ! folder should use NGUARD instead of these.
  integer, save :: gr_iguard = NGUARD
  integer, save :: gr_jguard = NGUARD
  integer, save :: gr_kguard = NGUARD

  integer, save :: gr_minRefine

  ! This is variable is for purely internal use with AMReX callbacks
  ! Do *not* use this directly
  logical, save :: gr_amrexDidRefinement = .FALSE.

!  integer,save,dimension(MDIM)::gr_bndOrder
!  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vartypes
  logical, save :: gr_justExchangedGC
!  logical, save :: gr_isolatedBoundaries
!  logical, save :: gr_useParticles,gr_refineOnParticleCount,gr_refineOnPdens
!  integer, save :: gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save :: gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs

!  logical, save :: gr_useEnergyDeposition
  real, save, dimension(LOW:HIGH,MDIM) :: gr_boxContainingLeafNodes

  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit
!  integer, save :: gr_oneRefLev=1 !! To be used with the multigrid
!  logical ,save :: gr_earlyBlockDistAdjustment
!  logical, save :: gr_monotone
!  integer, save :: gr_intpol

  integer, save :: gr_sanitizeDataMode
!  integer, save :: gr_sanitizeVerbosity

  integer, save :: gr_numRefineVars, gr_numRefineVarsMax
!  integer,allocatable,dimension(:) ,save :: gr_refineVars
  real,    save :: gr_globalDomain(LOW:HIGH, MDIM)
  integer, save :: gr_domainBC(2, MDIM)
!  integer,save,dimension(2,MDIM) :: gr_blkBC
  logical, save :: gr_dirIsAngular(MDIM)
!  logical, save :: gr_geometryOverride
  character(len=MAX_STRING_LENGTH) :: gr_str_geometry
  integer, save :: gr_refine_var(MAXREFVARS)
  real,    save :: gr_refine_cutoff(MAXREFVARS)
  real,    save :: gr_derefine_cutoff(MAXREFVARS)
  real,    save :: gr_refine_filter(MAXREFVARS)
  integer, save :: gr_globalNumBlocks !

  logical, save :: gr_doFluxCorrection

  logical, save :: gr_enableMaskedGCFill

  integer, save :: gr_interpolator

#ifdef GRID_WITH_MONOTONIC
  integer, save :: gr_intpolStencilWidth
#else
  ! The following was appropriate for Paramesh3f with native interpolation
  integer, parameter :: gr_intpolStencilWidth = 1
#endif
!
!#ifdef FLASH_PARTICLES
!  integer, save :: gr_maxParticlesPerProc
!#endif
!
!
!#ifdef FL_NON_PERMANENT_GUARDCELLS
!  integer,save :: gr_blkPtrRefCount, gr_blkPtrRefCount_fc
!  integer,save :: gr_lastBlkPtrGotten, gr_lastBlkPtrGotten_fc
!  logical,save,dimension(NUNK_VARS) :: gr_ccMask
!#if(NFACE_VARS>0)
!  logical,save,dimension(3,NFACE_VARS) :: gr_fcMask
!#else
!  logical :: gr_fcMask
!#endif
!
!#endif
!
  integer, save :: gr_lrefineDel
  logical, save :: gr_enforceMaxRefinement
!  logical,save :: gr_lrefineMaxRedDoByLogR, gr_lrefineMaxRedDoByTime
!  real,save :: gr_lrefineMaxRedRadiusSq
!  real,save :: gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
!  real,save :: gr_lrefineMaxRedTimeScale, gr_lrefineMaxRedTRef, gr_lrefineMaxRedLogBase
!  integer,save :: gr_restrictAllMethod
!
!  integer, save :: gr_lrefineMinInit
!  real, save, dimension(LOW:HIGH,MDIM) :: gr_region
!
!  logical,save :: gr_lrefinemaxByTime
!  real, save :: gr_lrefmaxTimes(GR_LREFMAXTIMES)
!  integer, save :: gr_lrefmaxTimeValues(GR_LREFMAXTIMES)
  logical, save ::  gr_gcellsUpToDate = .false.
!
!  logical, save :: gr_reduceGcellFills = .false.
!
!  real,allocatable,dimension(:) :: gr_error
end Module Grid_data
