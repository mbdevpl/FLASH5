!!****if* source/Grid/GridMain/AMR/Amrex/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!  Grid_init()
!!
!! DESCRIPTION
!!  Initialize the Grid unit and set initial values to all variables declared 
!!  in Grid_data
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS 
!!
!!  nblockx [INTEGER] 
!!     num initial blocks in x dir
!!  nblocky [INTEGER] 
!!     num initial blocks in y dir   
!!  nblockz [INTEGER] 
!!     num initial blocks in z dir   
!!  lrefine_max [INTEGER] 
!!      maximum AMR refinement level
!!  lrefine_min [INTEGER] 
!!      minimum AMR refinement level
!!  nrefs [INTEGER] 
!!      refine/derefine AMR grid every nrefs timesteps
!!
!!  refine_var_1 [INTEGER] 
!!     indicates first variable on which to refine
!!  refine_cutoff_1 [REAL] 
!!      threshold value to trigger refinement for refine_var_1
!!  derefine_cutoff_1 [REAL]
!!      threshold value to trigger derefinement for refine_var_1
!!  refine_filter_1 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_1
!!
!!  refine_var_2 [INTEGER] 
!!     indicates second variable on which to refine
!!  refine_cutoff_2 [REAL] 
!!      threshold value to trigger refinement for refine_var_2
!!  derefine_cutoff_2 [REAL]
!!      threshold value to trigger derefinement for refine_var_2
!!  refine_filter_2 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_2
!!
!!  refine_var_3 [INTEGER] 
!!     indicates third variable on which to refine (if needed)
!!  refine_cutoff_3 [REAL] 
!!      threshold value to trigger refinement for refine_var_3
!!  derefine_cutoff_3 [REAL]
!!      threshold value to trigger derefinement for refine_var_3
!!  refine_filter_3 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_3
!!
!!  refine_var_4 [INTEGER] 
!!     indicates fourth variable on which to refine (if needed)
!!  refine_cutoff_4 [REAL] 
!!      threshold value to trigger refinement for refine_var_4
!!  derefine_cutoff_4 [REAL]
!!      threshold value to trigger derefinement for refine_var_4
!!  refine_filter_4 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_4
!!
!!  flux_correct [BOOLEAN]
!!     turns on or off flux correction
!! small  [REAL]
!!   Generic small value that can be used as floor where needed
!! smlrho [REAL]  
!!   Cutoff value for density    
!! smallp [REAL]  
!!   Cutoff value for pressure
!! smalle [REAL]  
!!   Cutoff value for energy
!! smallt [REAL]  
!!   Cutoff value for temperature
!! smallu [REAL]  
!!   Cutoff value for velocity
!! smallx [REAL]  
!!   Cutoff value for abundances
!! eosMode[STRING]
!!   determines which variables to calculate from the ones
!!   defined. Possible values are "dens_ie", "dens_pres" and "dens_temp"
!! interpol_order [INTEGER]
!!   Order of interpolation, used in Paramesh2 "monotonic" interpolation
!!   for mesh prolongation
!! grid_monotone_hack [BOOLEAN]
!!   If .true., apply radical monotonicity constraints to interpolants,
!!   i.e., completely flatten them if they violate monotonicity.  Used
!!   in Paramesh2 "quadratic_cartesian" interpolation for mesh prolongation.
!! earlyBlockDistAdjustment [BOOLEAN]
!!   If .true., let Paramesh redistribute blocks
!!   across processors early, so that the block distribution chosen by
!!   Paramesh will be in effect when time evolution begins after restart.
!!   If earlyBlockDistAdjustment is .false., the block distribution enacted
!!   by the IO unit when it read a checkpoint file will normally still be
!!   in effect when time evolution begins after a restart.
!!   This flag is ignored if not restarting from a checkpoint.
!! 
!! lrefine_del [INTEGER]
!! gr_lrefineMaxRedDoByTime [BOOLEAN]
!! gr_lrefineMaxRedTRef [REAL]
!! gr_lrefineMaxRedTimeScale [REAL]
!! gr_lrefineMaxRedLogBase [REAL]
!!
!! gr_lrefineMaxRedDoByLogR [BOOLEAN]
!! gr_lrefineMaxRedRadiusFact [REAL]
!! x_refine_center [REAL]
!! y_refine_center [REAL]
!! z_refine_center [REAL]
!!
!! gr_restrictAllMethod [INTEGER]
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_init()
  use iso_c_binding,               ONLY : c_loc, c_null_ptr

  use amrex_base_module,           ONLY : amrex_spacedim
  use amrex_bc_types_module,       ONLY : amrex_bc_int_dir, &
                                          amrex_bc_ext_dir

  use Grid_data
  use Grid_interface,              ONLY : Grid_getDeltas, &
                                          Grid_getMaxRefinement
  use gr_amrexInterface,           ONLY : gr_amrexInit
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use Driver_interface,            ONLY : Driver_abortFlash, &
                                          Driver_getMype, &
                                          Driver_getNumProcs, &
                                          Driver_getComm
  use Logfile_interface,           ONLY : Logfile_stampMessage
  use Simulation_interface,        ONLY : Simulation_mapStrToInt, &
                                          Simulation_getVarnameType

  implicit none

  integer :: i, var

  character(len=MAX_STRING_LENGTH) :: paramString
  character(len=MAX_STRING_LENGTH) :: refVarname, refVarString
  character(len=MAX_STRING_LENGTH) :: refCutoffName, refCutOffString
  character(len=MAX_STRING_LENGTH) :: derefCutoffName, derefCutOffString
  character(len=MAX_STRING_LENGTH) :: refFiltername, refFilterString
  character(len=MAX_STRING_LENGTH) :: xl_bcString, xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString, yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString, zr_bcString
  character(len=MAX_STRING_LENGTH) :: eosModeString
  character(len=MAX_STRING_LENGTH) :: interpolatorString
  integer :: refVar
  integer :: nonrep

 
  call Driver_getMype(GLOBAL_COMM, gr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, gr_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM, gr_globalComm)

  call Driver_getMype(MESH_COMM, gr_meshMe)
  call Driver_getNumProcs(MESH_COMM, gr_meshNumProcs)
  call Driver_getComm(MESH_COMM, gr_meshComm)

  call Driver_getMype(MESH_ACROSS_COMM, gr_meshAcrossMe)
  call Driver_getNumProcs(MESH_ACROSS_COMM, gr_meshAcrossNumProcs)
  call Driver_getComm(MESH_ACROSS_COMM, gr_meshAcrossComm)

#ifdef GRID_WITH_MONOTONIC
  if (NGUARD < 4) then
     if (gr_meshMe==MASTER_PE) then
        print*,'Grid_init: Monotonic grid interpolation requires at least 4 layers of guard cells.'
        print*,' However, NGUARD is only ', NGUARD
        print*," Maybe you want to setup with '-gridinterpolation=native',"
        print*," or make sure that NGUARD is set correctly in Config file."
        call Driver_abortFlash("Please setup with '-gridinterpolation=native', or change NGUARD.")
     end if
  endif
#endif

!------------------------------------------------------------------------------
! Load into local Grid variables all runtime parameters needed by gr_initGeometry
!------------------------------------------------------------------------------
  call RuntimeParameters_get("geometry",gr_str_geometry)
  call RuntimeParameters_mapStrToInt(gr_str_geometry, gr_geometry)

  !get the boundary conditions stored as strings in the flash.par file
  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
  call RuntimeParameters_get("zr_boundary_type", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString, gr_domainBC(LOW, IAXIS))
  call RuntimeParameters_mapStrToInt(xr_bcString, gr_domainBC(HIGH,IAXIS))
  call RuntimeParameters_mapStrToInt(yl_bcString, gr_domainBC(LOW, JAXIS))
  call RuntimeParameters_mapStrToInt(yr_bcString, gr_domainBC(HIGH,JAXIS))
  call RuntimeParameters_mapStrToInt(zl_bcString, gr_domainBC(LOW, KAXIS))
  call RuntimeParameters_mapStrToInt(zr_bcString, gr_domainBC(HIGH,KAXIS))

!------------------------------------------------------------------------------
! FLASH inits geometry first as it can change runtime parameters
!------------------------------------------------------------------------------
  ! Determine the geometries of the individual dimensions, and scale
  ! angle value parameters that are expressed in degrees to radians.
  call gr_initGeometry()

!----------------------------------------------------------------------------------
! Initialize AMReX
!----------------------------------------------------------------------------------
  call gr_amrexInit()

  ! Check whether the dimensionality of the AMReX library matches the FLASH configuration
#ifdef __INTEL_COMPILER
#define SUBASSERT(asser) call abo(#asser)
#else
#define SUBASSERT(asser) call abo('asser')
#endif
#define ASSERT(assertion) if (.NOT.(assertion)) SUBASSERT(assertion)

  ASSERT(NDIM==amrex_spacedim)
  ASSERT(N_DIM==amrex_spacedim)


  ! Save BC information for AMReX callbacks
  lo_bc_amrex(:, :) = amrex_bc_int_dir
  hi_bc_amrex(:, :) = amrex_bc_int_dir
  do i = 1, NDIM
     select case(gr_domainBC(LOW, i))
     case(PERIODIC)
        lo_bc_amrex(i, :) = amrex_bc_int_dir
     case default
        lo_bc_amrex(i, :) = amrex_bc_ext_dir
     end select

     select case(gr_domainBC(HIGH, i))
     case(PERIODIC)
        hi_bc_amrex(i, :) = amrex_bc_int_dir
     case default
        hi_bc_amrex(i, :) = amrex_bc_ext_dir
     end select
  end do

!----------------------------------------------------------------------------------
! Store interface-accessible data as local Grid data variables for optimization
!----------------------------------------------------------------------------------
  !Store computational domain limits in a convenient array.  Used later in Grid_getBlkBC.
  call Grid_getDomainBoundBox(gr_globalDomain)

  call Grid_getMaxRefinement(gr_lRefineMax, mode=1)
  gr_maxRefine = gr_lRefineMax

  call Grid_getDeltas(gr_lRefineMax, gr_minCellSizes)
  gr_minCellSize = gr_minCellSizes(IAXIS)
#if NDIM >= 2
  if (.not. gr_dirIsAngular(JAXIS)) then
    gr_minCellSize = min(gr_minCellSize, gr_minCellSizes(JAXIS))
  end if
#if NDIM == 3
  if (.not. gr_dirIsAngular(KAXIS)) then
     gr_minCellSize = min(gr_minCellSize, gr_minCellSizes(KAXIS))
  end if
#endif
#endif

  gr_allPeriodic = .TRUE.
  do i = 1, NDIM
     if(gr_domainBC(LOW,  i) /= PERIODIC)     gr_allPeriodic = .FALSE.
     if(gr_domainBC(HIGH, i) /= PERIODIC)     gr_allPeriodic = .FALSE.
  end do

  ! DEV: TODO: Is there a FLASH Grid_*/gr_* interface function to get this?  If
  ! so, set it up and have it get the value from AMReX
  call RuntimeParameters_get("nrefs", gr_nrefs)

!----------------------------------------------------------------------------------
! Setup all remaining local Grid data variables
!----------------------------------------------------------------------------------
!  call RuntimeParameters_get("geometryOverride",gr_geometryOverride)

  call RuntimeParameters_get("convertToConsvdForMeshCalls", &
                             gr_convertToConsvdForMeshCalls)
  call RuntimeParameters_get("convertToConsvdInMeshInterp", &
                             gr_convertToConsvdInMeshInterp)
  if (gr_convertToConsvdForMeshCalls) then
    if (gr_meshMe == MASTER_PE) then
      print*,'WARNING: convertToConsvdForMeshCalls not implemented for AMReX.'
      print*,'         Using convertToConsvdInMeshInterp instead.'
    end if
    call Logfile_stampMessage( &
           'WARNING: convertToConsvdForMeshCalls not implemented for AMReX.' // &
           '         Using convertToConsvdInMeshInterp instead.')

    gr_convertToConsvdInMeshInterp = .TRUE.
    gr_convertToConsvdForMeshCalls = .FALSE.
  end if

  call RuntimeParameters_get("enableMaskedGCFill", gr_enableMaskedGCFill)
  
  ! DEV: FIXME Temporarily fix so that sims halt on error.  Allow user
  ! specified value once this functionality is fully implemented
  gr_sanitizeDataMode = 4
!  call RuntimeParameters_get("gr_sanitizeDataMode",  gr_sanitizeDataMode)
!  call RuntimeParameters_get("gr_sanitizeVerbosity", gr_sanitizeVerbosity)

  call RuntimeParameters_get("smalle", gr_smalle)
  call RuntimeParameters_get("smlrho", gr_smallrho)
!  call RuntimeParameters_get("smallx",gr_smallx) !

!!  call RuntimeParameters_get("grid_monotone_hack", gr_monotone) ! for "quadratic_cartesian" interpolation
!  call RuntimeParameters_get("interpol_order",gr_intpol) ! for "monotonic" interpolation
  call RuntimeParameters_get("amrexInterpolator", interpolatorString)
  call RuntimeParameters_mapStrToInt(interpolatorString, gr_interpolator)
  if (gr_interpolator == NONEXISTENT) then
    call Driver_abortFlash("[Grid_init] Unknown amrexInterpolator runtime parameter value")
  end if

#ifdef GRID_WITH_MONOTONIC
  gr_intpolStencilWidth = 2     !Could possibly be less if gr_intpol < 2  - KW
#endif

!  call RuntimeParameters_get("bndPriorityOne",gr_bndOrder(1))
!  call RuntimeParameters_get("bndPriorityTwo",gr_bndOrder(2))
!  call RuntimeParameters_get("bndPriorityThree",gr_bndOrder(3))

!  call RuntimeParameters_get("refine_on_particle_count",gr_refineOnParticleCount)

!  call RuntimeParameters_get("min_particles_per_blk",gr_minParticlesPerBlk)
!  call RuntimeParameters_get("max_particles_per_blk",gr_maxParticlesPerBlk)

  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosMode)

  call RuntimeParameters_get("eosModeInit", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosModeInit)

!  call RuntimeParameters_get("earlyBlockDistAdjustment", gr_earlyBlockDistAdjustment)
  gr_justExchangedGC = .FALSE.

  call gr_initSpecific()

  !! This section of the code identifies the variables to used in
  !! the refinement criterion. If a variable is a refinement variable
  !! then the corresponding refinement/derefinement cutoff and filter
  !! values also have to be fetched. The config file defines
  !! refinement variables as strings, names as "refine_var_1",
  !! "refine_var_2" etc, with the current maximum being 4. The
  !! general utility routine takes the base "refine_var_" and appends
  !! the index at the end of the string to generate the parameter
  !! name and the routine Simulation_mapStrToInt finds its index into UNK.

  call RuntimeParameters_get("refine_var_count", gr_numRefineVarsMax)
  gr_refine_var = NONEXISTENT
  gr_numRefineVars = 0

  refVarName = 'refine_var_'
  refCutoffName = 'refine_cutoff_'
  derefCutoffName ='derefine_cutoff_'
  refFilterName = 'refine_filter_'

  do i = 1, gr_numRefineVarsMax
     call concatStringWithInt(refVarName, i, refVarString)
     call RuntimeParameters_get(refVarString, paramString)
     if(paramString /= "none") then
        do ! not a real loop
           call Simulation_mapStrToInt(paramString, refVar, MAPBLOCK_UNK)
           if(refVar <= 0)  EXIT
           nonrep = 0
           ! DEV: FIXME This is segfaulting, but isn't needed right now
!           call Grid_getVarNonRep(MAPBLOCK_UNK, refVar, nonrep)
           if(nonrep > 0) then
             refVar = 0
             exit
           end if

           gr_numRefineVars = gr_numRefineVars + 1
           gr_refine_var(gr_numRefineVars) = refVar
           call concatStringWithInt(refCutoffName, gr_numRefineVars, &
                                    refCutoffString)
           call concatStringWithInt(derefCutoffName, gr_numRefineVars, &
                                    derefCutOffString)
           call concatStringWithInt(refFilterName, gr_numRefineVars, &
                                    refFilterString)
           call RuntimeParameters_get(refCutoffString, &
                                      gr_refine_cutoff(gr_numRefineVars))
           call RuntimeParameters_get(derefCutoffString, &
                                      gr_derefine_cutoff(gr_numRefineVars))
           call RuntimeParameters_get(refFilterString, &
                                      gr_refine_filter(gr_numRefineVars))
           exit ! told you it wasnt a real loop
        end do
        if(refVar <= 0) then
           if(gr_globalMe == MASTER_PE) &
              print*, 'WARNING: Unrecognized or non-replicated variable name in refine_var_',i,' treating it as "none"'
           call Logfile_stampMessage( &
              'WARNING: Unrecognized or non-replicatedvariable name in refine_var, treating it as "none"')
        end if
     end if
  end do

  if (gr_numRefineVars == 0) then
     if (gr_meshMe == MASTER_PE) then
        print*,'WARNING : Adaptive Grid did not find any refinement variables'
     end if
     call Logfile_stampMessage("WARNING : Adaptive Grid did not find any variable to refine")
  end if

  gr_enforceMaxRefinement = .FALSE.

  call RuntimeParameters_get('lrefine_min', gr_minRefine)
  call RuntimeParameters_get("lrefine_del", gr_lrefineDel)

!  call RuntimeParameters_get("gr_lrefineMaxRedDoByLogR", gr_lrefineMaxRedDoByLogR)
!  call RuntimeParameters_get("gr_lrefineMaxRedRadiusFact", gr_lrefineMaxRedRadiusSq)
!  gr_lrefineMaxRedRadiusSq = gr_lrefineMaxRedRadiusSq * gr_lrefineMaxRedRadiusSq
!  call RuntimeParameters_get("x_refine_center", gr_lrefineCenterI)
!  call RuntimeParameters_get("y_refine_center", gr_lrefineCenterJ)
!  call RuntimeParameters_get("z_refine_center", gr_lrefineCenterK)
!
!  call RuntimeParameters_get("gr_lrefineMaxRedDoByTime", gr_lrefineMaxRedDoByTime)
!  if (gr_lrefineMaxRedDoByTime) gr_enforceMaxRefinement = .TRUE.
!  call RuntimeParameters_get("gr_lrefineMaxRedTimeScale", gr_lrefineMaxRedTimeScale)
!  call RuntimeParameters_get("gr_lrefineMaxRedLogBase", gr_lrefineMaxRedLogBase)
!  call RuntimeParameters_get("gr_lrefineMaxRedTRef", gr_lrefineMaxRedTRef)
!
!  call RuntimeParameters_get("gr_lrefineMaxByTime", gr_lrefineMaxByTime)
!  if (gr_lrefineMaxByTime) then
!     gr_enforceMaxRefinement = .TRUE.
!
!     do i = 1, GR_LREFMAXTIMES
!        refVarName = "gr_lrefmaxTime_"
!        call concatStringWithInt(refVarName,i,refVarString)
!        call RuntimeParameters_get(refVarString, gr_lrefmaxTimes(i))
!
!        refVarName = "gr_lrefmaxTimeValue_"
!        call concatStringWithInt(refVarName,i,refVarString)
!        call RuntimeParameters_get(refVarString, gr_lrefmaxTimeValues(i))
!
!        if(i > 1) then 
!           if(gr_lrefmaxTimes(i) > 0.0 .and. &
!                gr_lrefmaxTimes(i) < gr_lrefmaxTimes(i-1)) then
!              if(gr_meshMe == MASTER_PE) then
!                 call Driver_abortFlash('[Grid_init] Custom lrefine_max times must be in increasing order')
!              end if
!           end if
!        end if
!
!        if(gr_lrefmaxTimes(i) > 0.0 .and. gr_lrefmaxTimeValues(i) < 0) then
!           call Driver_abortFlash('[Grid_init] If custom lrefine_max time is set, the value must also be set')
!        end if
!     end do
!  end if

!  call RuntimeParameters_get("gr_restrictAllMethod", gr_restrictAllMethod)

    call RuntimeParameters_get("flux_correct", gr_doFluxCorrection)

!#ifdef FLASH_PARTICLES
!  call RuntimeParameters_get('useParticles',gr_useParticles)
!  call RuntimeParameters_get('pt_maxPerProc',gr_maxParticlesPerProc)
!#else
!  gr_useParticles=.false.
!#endif
!
!#ifdef FLASH_EDEP
!  call RuntimeParameters_get('useEnergyDeposition', gr_useEnergyDeposition)
!  gr_useParticles = gr_useEnergyDeposition
!#else
!  gr_useEnergyDeposition = .false.
!#endif
!
!#ifdef FLASH_GRID_PARTICLES
!  call RuntimeParameters_get('useProtonImaging',useProtonImaging)
!  if (useProtonImaging) then
!      gr_useParticles=.true.
!  end if
!  call RuntimeParameters_get('useProtonEmission',useProtonEmission)
!  if (useProtonEmission) then
!      gr_useParticles=.true.
!  end if
!#endif

  !Check if there are gravitational isolated boundary conditions
  !in order to determine which solvers to intialize.
!  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
!  gr_isolatedBoundaries = (grav_boundary_type=="isolated")
!
  do i = UNK_VARS_BEGIN,UNK_VARS_END
!     gr_vars(i)=i
     call Simulation_getVarnameType(i, gr_vartypes(i))
  end do

!#ifdef FL_NON_PERMANENT_GUARDCELLS
!  gr_blkPtrRefCount = 0 
!  gr_blkPtrRefCount_fc = 0
!  gr_lastBlkPtrGotten = 0 
!  gr_lastBlkPtrGotten_fc = 0
!#endif
!  call gr_setDataStructInfo()
!  call gr_bcInit()

  ! DEVNOTE: Is this Paramesh-specific.  Defined in paramesh
  !Initialize grid arrays used by IO
!  allocate(gr_nToLeft(0:gr_meshNumProcs-1))
!  allocate(gr_gid(nfaces+nchild+1, MAXBLOCKS))

  !Only call the particle initialization routines when
  !we are using particles.
!  if (gr_useParticles .eqv. .true. ) then
!     call gr_ptInit()
!     call gr_ptMapInit()
!  endif

!  call gr_sbInit()

    ! Reduce guard cell fills
!  call RuntimeParameters_get ("reduceGcellFills", gr_reduceGcellFills)

!  gr_region=0.0
  
! DEVNOTE: What to do with these?
!#ifndef BSS_GRID_ARRAYS
!# if NSCRATCH_GRID_VARS > 0
!  allocate(scratch(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
!       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc+1,&
!       gr_kLoGc:gr_kHiGc+1,MAXBLOCKS))
!# else
!  allocate(scratch(1,1,1,1,1))
!# endif
!
!# if NSCRATCH_CENTER_VARS > 0
!  allocate(scratch_ctr(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
!       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
!       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
!# else
!  allocate(scratch_ctr(1,1,1,1,1))
!# endif
!
!# if(NSCRATCH_FACEX_VARS>0)  
!  allocate(scratch_facevarx( SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
!       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
!       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
!# else
!  allocate(scratch_facevarx(1,1,1,1,1))
!# endif
!
!# if(NSCRATCH_FACEY_VARS>0)  
!  allocate(scratch_facevary( SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
!       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
!       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
!# else
!  allocate(scratch_facevary(1,1,1,1,1))
!# endif  
!
!# if(NSCRATCH_FACEZ_VARS>0)
!  allocate(scratch_facevarz( SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
!       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
!       gr_kLoGc:gr_kHiGc+K3D,MAXBLOCKS) )
!# else
!  allocate(scratch_facevarz(1,1,1,1,1))
!# endif

!  allocate(gr_xflx(NFLUXES,2,NYB,NZB,MAXBLOCKS))
!  allocate(gr_yflx(NFLUXES,NXB,2,NZB,MAXBLOCKS))
!  allocate(gr_zflx(NFLUXES,NXB,NYB,2,MAXBLOCKS))
!  
!# ifdef FLASH_HYDRO_UNSPLIT
!#  if NDIM >= 2
!  allocate(gr_xflx_yface(NFLUXES,2:NXB, 2   ,NZB  ,MAXBLOCKS))
!  allocate(gr_yflx_xface(NFLUXES,2    ,2:NYB,NZB  ,MAXBLOCKS))
!#   if NDIM == 3
!  allocate(gr_xflx_zface(NFLUXES,2:NXB,NYB  , 2   ,MAXBLOCKS))
!  allocate(gr_yflx_zface(NFLUXES,NXB,  2:NYB, 2   ,MAXBLOCKS))
!  allocate(gr_zflx_xface(NFLUXES, 2 ,NYB    ,2:NZB,MAXBLOCKS))
!  allocate(gr_zflx_yface(NFLUXES,NXB, 2     ,2:NZB,MAXBLOCKS))
!#   endif
!#  endif
!# endif

!#endif

  if(gr_meshMe == MASTER_PE) call printRefinementInfo()

contains

  subroutine printRefinementInfo()
    
    implicit none

    integer           :: ilev = 0
    integer           :: n = 0
    real              :: del(MDIM) = 0.0d0
    character(len=20) :: fmtStr = ""
    character(len=2)  :: colHdr(MDIM) = (/'dx', 'dy', 'dz'/)

    write(*,*) 'Grid_init: resolution based on runtime params:'
    write(*,'(A9,3(A12:4x))')  'lrefine', (colHdr(n),n=1,NDIM)
    do ilev = 1, gr_lRefineMax
       call Grid_getDeltas(ilev, del)
       if (maxval(del(IAXIS:NDIM)) .GT. 999999999999.999) then
          fmtStr = '(I7,2x,1P,3G16.3)'
       else if (minval(del(IAXIS:NDIM)) .LE. 0.0009) then
          fmtStr = '(I7,2x,1P,3G16.3)'
       else
          fmtStr = '(I7,2x,3F16.3)'
       end if

       write(*,fmtStr) ilev, (del(n),n=1,NDIM)
    end do
  end subroutine printRefinementInfo

  subroutine abo(msg)
    character(len=*) msg
    call Driver_abortFlash("Failed assertion in Grid_init: "//msg)
  end subroutine abo
end subroutine Grid_init
