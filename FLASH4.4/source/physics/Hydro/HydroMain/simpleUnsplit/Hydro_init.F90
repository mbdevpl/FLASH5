!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/HLL/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  Hydro_init()
!!  
!!
!! DESCRIPTION
!! 
!!  This routine initializes unit scope variables which are typically the runtime parameters.
!!  The routine must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!
!!***

Subroutine Hydro_init()

  use Hydro_data
  use Driver_interface,            ONLY : Driver_abortFlash, Driver_getMype, &
                                          Driver_getComm,                    &
                                          Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stampMessage, &
                                          Logfile_stampVarMask, &
                                          Logfile_stamp
  use Grid_interface,              ONLY : Grid_setFluxHandling


  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer :: i
  logical :: threadBlockListBuild, threadWithinBlockBuild

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)
  call Driver_getComm(MESH_COMM,hy_meshComm)

  call RuntimeParameters_get("useHydro",hy_useHydro)

  call RuntimeParameters_get("cfl", hy_cfl)
  hy_cfl_original = hy_cfl

  call RuntimeParameters_get("UnitSystem",          hy_units)
  call RuntimeParameters_get("order",               hy_order)
  call RuntimeParameters_get("hy_3Torder",          hy_3Torder)
  if (hy_3Torder == -1) hy_3Torder = hy_order
  call RuntimeParameters_get("transOrder",          hy_transOrder)
  if ((NGUARD <= 4) .and. (hy_order > 3)) then
     call Driver_abortFlash&
          ("[Hydro_init]: Hydro requires more guardcells for the given hy_order method.")
  endif
  hy_useVaryingCFL = .false.
  call RuntimeParameters_get("use_hybridOrder",     hy_useHybridOrder)
  if (hy_useHybridOrder) hy_useVaryingCFL = .true.
  call RuntimeParameters_get("hybridOrderKappa",    hy_hybridOrderKappa)

#ifdef BDRY_VAR
  if (.not. hy_useVaryingCFL) hy_useVaryingCFL = .true.
#endif

#ifdef FLASH_GRID_PARAMESH
  if ((NGUARD > 4) .and. (NXB < 2*NGUARD)) then
     call Driver_abortFlash&
          ("[Hydro_init]: Hydro requires larger NXB, etc. for the given number of guardcells.")
  endif
#endif


  call RuntimeParameters_get("entropy",             hy_entropy)
  call RuntimeParameters_get('entropyFixMethod',    hy_entropyFixMethod_str)
  call RuntimeParameters_get("charLimiting",        hy_charLimiting)
  call RuntimeParameters_get("eintSwitch",          hy_eswitch)
  call RuntimeParameters_get("slopeLimiter",        hy_limiter_str)
  call RuntimeParameters_get("smlrho",              hy_smalldens)
  call RuntimeParameters_get("smallp",              hy_smallpres)
  call RuntimeParameters_get("smallu",              hy_smallu)
  call RuntimeParameters_get('irenorm',             hy_irenorm)
  call RuntimeParameters_get('RiemannSolver',       hy_RiemannSolver_str)
  call RuntimeParameters_get('shockDetect',         hy_shockDetectOn)
  call RuntimeParameters_get("use_steepening",      hy_ContactSteepening)
  call RuntimeParameters_get("EOSforRiemann",       hy_EOSforRiemann)
  call RuntimeParameters_get("use_flattening",      hy_flattening)
  call RuntimeParameters_get("use_upwindTVD",       hy_upwindTVD)
  if (NGUARD <= 4) hy_upwindTVD = .false.
  call RuntimeParameters_get("use_avisc",           hy_use_avisc)
  call RuntimeParameters_get("cvisc",               hy_cvisc)
  call RuntimeParameters_get("use_GravPotUpdate",   hy_useGravPotUpdate)
  call RuntimeParameters_get("updateHydroFluxes",   hy_updateHydrofluxes)
  call RuntimeParameters_get("addThermalFlux",      hy_addThermalFlux)
  call RuntimeParameters_get("conserveAngMom",      hy_conserveAngMom)  
  if (NDIM == 3) then
     call RuntimeParameters_get("use_3dFullCTU",    hy_use3dFullCTU)
  else
     hy_use3dFullCTU = .false.
  endif


  !! Gravity -------------------------------------------------------------------
  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
  if (hy_useGravity) then
     call RuntimeParameters_get("use_gravHalfUpdate", hy_useGravHalfUpdate)
     if (hy_useGravHalfUpdate) then
        call RuntimeParameters_get("use_gravConsv",   hy_gravConsv)
     else
        hy_gravConsv = .false.
     endif

     if (hy_useGravPotUpdate .and. hy_gravConsv) then
        call Driver_abortFlash&
          ("[Hydro_init]: Gravity update using Poisson solver only support primitive update."//&
           "'Please use hy_gravConsv = .false.'")
     endif
  endif
#endif


  !! Non-ideal diffusions ------------------------------------------------------
  call RuntimeParameters_get("useViscosity",    hy_useViscosity)
  call RuntimeParameters_get("useConductivity", hy_useConductivity)

  hy_useDiffuse = .false.
  if (hy_useViscosity .or. hy_useConductivity) then
     hy_useDiffuse = .true.
  endif


  !! Entropy Fix ---------------------------------------------------------------
#if(0)
  if (hy_entropy) then
     if(trim(hy_entropyFixMethod_str) == "HartenHyman" .or. &
        trim(hy_limiter_str) == "hartenhyman") then
        hy_entropyFixMethod = HARTENHYMAN
     elseif(trim(hy_entropyFixMethod_str) == "Harten" .or. &
            trim(hy_limiter_str) == "harten") then
        hy_entropyFixMethod = HARTEN
     endif
  endif
#endif


  !! Slope Limiter -------------------------------------------------------------
#if(0)
  if(trim(hy_limiter_str) == "minmod" .or. &
     trim(hy_limiter_str) == "MINMOD") then
     hy_limiter = MINMOD
  else if(trim(hy_limiter_str) == "mc" .or. &
          trim(hy_limiter_str) == "MC" ) then
     hy_limiter = MC
  else if (trim(hy_limiter_str) == "hybrid" .or. &
           trim(hy_limiter_str) == "HYBRID" ) then
     hy_limiter = HYBRID
  else if (trim(hy_limiter_str) == "vanLeer" .or. &
           trim(hy_limiter_str) == "VANLEER") then
     hy_limiter = VANLEER
  else if (trim(hy_limiter_str) == "limited" .or. &
           trim(hy_limiter_str) == "LIMITED") then
     hy_limiter = LIMITED
     call RuntimeParameters_get('LimitedSlopeBeta', hy_LimitedSlopeBeta)
  else
     call Driver_abortFlash&
          ("[Hydro_init]: The hy_limter of unknown type! It should be one of" // &
           "'minmod','mc', 'vanLeer', 'hybrid' or 'limited'.")
  end if
#endif


  !! Riemann Solver ------------------------------------------------------------
  hy_hybridRiemannOnly = .false.

  if(trim(hy_RiemannSolver_str) == "Roe" .or. &
     trim(hy_RiemannSolver_str) == "roe" .or. &
     trim(hy_RiemannSolver_str) == "ROE" ) then
     hy_RiemannSolver = ROE
  elseif (trim(hy_RiemannSolver_str) == "hll" .or. &
          trim(hy_RiemannSolver_str) == "HLL" ) then
     hy_RiemannSolver = HLL
#if(0)
  elseif (trim(hy_RiemannSolver_str) == "hllc" .or. &
          trim(hy_RiemannSolver_str) == "HLLC" ) then
     hy_RiemannSolver = HLLC
  elseif (trim(hy_RiemannSolver_str) == "marquina" .or. &
          trim(hy_RiemannSolver_str) == "Marquina" ) then
     hy_RiemannSolver = MARQ
#endif
  elseif (trim(hy_RiemannSolver_str) == "LocalLaxFriedrichs" .or. &
          trim(hy_RiemannSolver_str) == "llf"  .or. &
          trim(hy_RiemannSolver_str) == "LLF" ) then
     hy_RiemannSolver = LLF
#if(0)
  elseif (trim(hy_RiemannSolver_str) == "HYBRID" .or. &
          trim(hy_RiemannSolver_str) == "hybrid"  .or. &
          trim(hy_RiemannSolver_str) == "Hybrid" ) then
     hy_RiemannSolver = HYBR
     if (.not.hy_shockDetectOn) then
        hy_hybridRiemannOnly = .true.
        hy_shockDetectOn = .true.
     endif
#endif
  else
     call Driver_abortFlash&
          ("[Hydro_init]: The Riemann Solver is of unknown type: " // &
           "Options are HLL or LLF.")
  endif

  !! Geometry ------------------------------------------------------------------
  call RuntimeParameters_get("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, hy_geometry)
  if (hy_geometry .NE. CARTESIAN )  then
     if (hy_meshME == MASTER_PE) print *, "[Hydro_init]: Using non-Cartesian Geometry!"
  endif
  call RuntimeParameters_get("flux_correct", hy_fluxCorrect)
  if (NDIM > 1) then
     if (hy_fluxCorrect) then
        if (hy_geometry == CARTESIAN) then
           call Grid_setFluxHandling('consv_flux_densities')
        else
           call Grid_setFluxHandling('consv_fluxes')
        endif
     end if
  end if

  call PhysicalConstants_get("ideal gas constant", hy_Rconst)

  !! For correct flux correction in non-Cartesian geometry----------------------
  do i = 1, NFLUXES
     hy_fluxCorVars(i) = i
  enddo
     

  !! System units --------------------------------------------------------------
  hy_xref = 1.0
  hy_vref = 1.0
  hy_dref = 1.0

  hy_mref = hy_xref*hy_vref
  hy_tref = hy_xref/hy_vref
  hy_eref = hy_vref*hy_vref
  hy_nref = hy_dref*hy_vref*hy_xref
  hy_pref = hy_dref*hy_vref*hy_vref
  hy_gref = hy_vref*hy_vref/hy_xref

  hy_qref = hy_vref*hy_vref/hy_Rconst
  hy_kref = hy_dref*hy_vref*hy_xref*hy_Rconst


  !! Allow selective guardcell fill calls ---------------------------------------
  hy_gcMaskSize = NUNK_VARS
  hy_gcMask = .TRUE.

#ifdef DFCF_VAR
  hy_gcMask(DFCF_VAR) = .FALSE.
#endif
#ifdef FLLM_VAR
  hy_gcMask(FLLM_VAR) = .FALSE.
#endif
#ifdef PIPE_VAR
  hy_gcMask(PIPE_VAR) = .FALSE.
#endif
#ifdef TITE_VAR
  hy_gcMask(TITE_VAR) = .FALSE.
#endif
#ifdef DBGS_VAR
  hy_gcMask(DBGS_VAR) = .FALSE.
#endif
#ifdef DENS_VAR
!!  hy_gcMask(DENS_VAR) = .FALSE.
#endif
#ifdef EELE_VAR
!!  hy_gcMask(EELE_VAR) = .FALSE.
#endif
#ifdef EINT_VAR
!!  hy_gcMask(EINT_VAR) = .FALSE.
#endif
#ifdef EION_VAR
!!  hy_gcMask(EION_VAR) = .FALSE.
#endif
#ifdef ENER_VAR
!!  hy_gcMask(ENER_VAR) = .FALSE.
#endif
#ifdef ERAD_VAR
!!  hy_gcMask(ERAD_VAR) = .FALSE.
#endif
#ifdef GAMC_VAR
  hy_gcMask(GAMC_VAR) = .TRUE.
#endif
#ifdef GAME_VAR
  hy_gcMask(GAME_VAR) = .FALSE.
#endif
#ifdef PELE_VAR
  hy_gcMask(PELE_VAR) = .FALSE.
#endif
#ifdef PION_VAR
  hy_gcMask(PION_VAR) = .FALSE.
#endif
#ifdef PRAD_VAR
  hy_gcMask(PRAD_VAR) = .FALSE.
#endif
#ifdef PRES_VAR
!!$  hy_gcMask(PRES_VAR) = .FALSE.
#endif
#ifdef TELE_VAR
  hy_gcMask(TELE_VAR) = .FALSE.
#endif
#ifdef TEMP_VAR
  hy_gcMask(TEMP_VAR) = .FALSE.
#endif
#ifdef TION_VAR
  hy_gcMask(TION_VAR) = .FALSE.
#endif
#ifdef TRAD_VAR
  hy_gcMask(TRAD_VAR) = .FALSE.
#endif
#ifdef VELX_VAR
!!  hy_gcMask(VELX_VAR) = .FALSE.
#endif
#ifdef VELY_VAR
!!  hy_gcMask(VELY_VAR) = .FALSE.
#endif
#ifdef VELZ_VAR
!!  hy_gcMask(VELZ_VAR) = .FALSE.
#endif
#ifdef VOLX_VAR
  hy_gcMask(VOLX_VAR) = .FALSE.
#endif
#ifdef VOLY_VAR
  hy_gcMask(VOLY_VAR) = .FALSE.
#endif
#ifdef VOLZ_VAR
  hy_gcMask(VOLZ_VAR) = .FALSE.
#endif

#if NSPECIES == 1
#ifdef SPECIES_BEGIN
  hy_gcMask(SPECIES_BEGIN) = .FALSE.
#endif
#endif

  call Logfile_stampVarMask(hy_gcMask, .FALSE., '[Hydro_init]', 'gcNeed')


  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", hy_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", hy_threadWithinBlock)

  if (hy_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadBlockList = .false.
  end if
  if (hy_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadWithinBlock = .false.
  end if


End Subroutine Hydro_init
