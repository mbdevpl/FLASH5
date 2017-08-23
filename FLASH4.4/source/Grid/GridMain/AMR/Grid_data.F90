!!****if* source/Grid/GridMain/paramesh/Grid_data
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
!!  This includes the global integer identifier for a block, the grid geometry information
!!  
!!  
!! 
!! CREATE AD:04/12/04
!!
!! 
!!   Defining data structures for storing paramesh related infomation.  
!!   including function for updating the grid information
!!
!! MODIFIED AD:05/19/04
!!   
!!***

Module Grid_data

  implicit none

#include "constants.h"
#include "Flash.h"

  real,save, allocatable, dimension(:,:) :: gr_delta
  integer, save :: gr_iguard = NGUARD
  integer, save :: gr_jguard = NGUARD 
  integer, save :: gr_kguard = NGUARD

  integer,save,dimension(MDIM)::gr_bndOrder
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vartypes
  logical, save :: gr_anyVarToConvert
  logical, save :: gr_justExchangedGC, gr_allPeriodic, gr_isolatedBoundaries
  logical, save :: gr_useParticles,gr_refineOnParticleCount,gr_refineOnPdens
  integer, save :: gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save :: gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs

  logical, save :: gr_useEnergyDeposition
  real, save, dimension(LOW:HIGH,MDIM) :: gr_boxContainingLeafNodes


  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  integer, save :: gr_oneRefLev=1 !! To be used with the multigrid
  integer ,save :: gr_nrefs
  logical ,save :: gr_convertToConsvdForMeshCalls
  logical ,save :: gr_convertToConsvdInMeshInterp
  logical ,save :: gr_earlyBlockDistAdjustment
  logical, save :: gr_monotone
  integer, save :: gr_intpol
  real, save :: gr_smallx

  integer, save :: gr_numRefineVars, gr_numRefineVarsMax
  integer,allocatable,dimension(:) ,save :: gr_refineVars
  real ,save :: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain
  integer,save,dimension(2,MDIM) :: gr_domainBC, gr_blkBC
  integer ,save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  character(len=MAX_STRING_LENGTH) :: gr_str_geometry
  integer ,save, dimension(MAXREFVARS) :: gr_refine_var    
  real,dimension(MAXREFVARS), save::gr_refine_cutoff,&
       gr_derefine_cutoff,gr_refine_filter
  real, save :: gr_smalle,gr_smallrho
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  integer, save :: gr_globalNumBlocks !


#ifdef GRID_WITH_MONOTONIC
  integer, save :: gr_intpolStencilWidth
#else
  ! The following is for Paramesh3f with native interpolation
  integer, parameter :: gr_intpolStencilWidth = 1
#endif

#ifdef FLASH_PARTICLES
  integer, save :: gr_maxParticlesPerProc
#endif

  integer, save :: gr_numDataStruct
  integer, save, dimension(NDATATYPES) :: gr_gridDataStruct,gr_gridDataStructSize

#ifdef FL_NON_PERMANENT_GUARDCELLS
  integer,save :: gr_blkPtrRefCount, gr_blkPtrRefCount_fc
  integer,save :: gr_lastBlkPtrGotten, gr_lastBlkPtrGotten_fc
  logical,save,dimension(NUNK_VARS) :: gr_ccMask
#if(NFACE_VARS>0)
  logical,save,dimension(3,NFACE_VARS) :: gr_fcMask
#else
  logical :: gr_fcMask
#endif

#endif

  integer,save :: gr_lrefineDel, gr_maxRefine, gr_minRefine
  logical,save :: gr_enforceMaxRefinement
  logical,save :: gr_lrefineMaxRedDoByLogR, gr_lrefineMaxRedDoByTime
  real,save :: gr_lrefineMaxRedRadiusSq
  real,save :: gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
  real,save :: gr_lrefineMaxRedTimeScale, gr_lrefineMaxRedTRef, gr_lrefineMaxRedLogBase
  integer,save :: gr_restrictAllMethod

  integer, save :: gr_lrefineMinInit
  real, save, dimension(LOW:HIGH,MDIM) :: gr_region

  logical,save :: gr_lrefinemaxByTime
  real, save :: gr_lrefmaxTimes(GR_LREFMAXTIMES)
  integer, save :: gr_lrefmaxTimeValues(GR_LREFMAXTIMES)
  logical, save ::  gr_gcellsUpToDate = .false.

  logical, save :: gr_reduceGcellFills = .false.

  logical, save :: gr_bcEnableApplyMixedGds
  real,allocatable,dimension(:) :: gr_error
end Module Grid_data
