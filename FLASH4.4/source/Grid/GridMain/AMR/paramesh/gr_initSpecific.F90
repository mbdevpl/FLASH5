!!****if* source/Grid/GridMain/paramesh/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!
!!  Grid_init()
!!           
!!
!! DESCRIPTION
!!  Initialize Grid_data
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

!!REORDER(5):scratch, scratch_ctr, scratch_facevar[xyz], gr_[xyz]flx
!!REORDER(5):gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
#include "constants.h"
subroutine gr_initSpecific()
  use Grid_data
  use gr_specificData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use paramesh_comm_data, ONLY : amr_mpi_meshComm
  use tree, ONLY : lrefine_min, lrefine_max, nfaces, nchild

  implicit none
  
  real :: dx, dy, dz
  real, dimension(NDIM) :: rnb
  integer :: i

  !get the initial grid layout
  call RuntimeParameters_get("nblockx", gr_nBlockX) !number of initial blks in x dir
  call RuntimeParameters_get("nblocky", gr_nBlockY) !number of initial blks in y dir
  call RuntimeParameters_get("nblockz", gr_nblockZ) !number of initial blks in z dir  

  call RuntimeParameters_get("enableMaskedGCFill", gr_enableMaskedGCFill)
  call RuntimeParameters_get("gr_sanitizeDataMode",  gr_sanitizeDataMode)
  call RuntimeParameters_get("gr_sanitizeVerbosity", gr_sanitizeVerbosity)

  amr_mpi_meshComm=gr_meshComm
  call Paramesh_init()
  call gr_initGeometry()

  !! calculating deltas for each level of 
  !! refinement and putting them in the
  !! delta variable
  dx = gr_imax - gr_imin
  dy = gr_jmax - gr_jmin
  dz = gr_kmax - gr_kmin
  rnb = 0.0
  rnb(1) = dx/(1.0*NXB*gr_nBlockX)
#if NDIM > 1
  rnb(2) = dy/(1.0*NYB*gr_nBlockY)
#endif
#if NDIM > 2
  rnb(3) = dz/(1.0*NZB*gr_nBlockZ)
#endif  
  do i = 1,gr_maxRefine
     gr_delta(1:NDIM,i) = rnb
     gr_delta(NDIM+1:,i) = 0.0
     rnb = rnb/2.0
  end do

  gr_minCellSizes(IAXIS) = (gr_imax - gr_imin) / &
       (gr_nblockX*NXB*2**(gr_maxRefine-1))
  gr_minCellSize = gr_minCellSizes(IAXIS)


  if (NDIM >= 2) then
     gr_minCellSizes(JAXIS) = (gr_jmax - gr_jmin) / &
          (gr_nblockY*NYB*2**(gr_maxRefine-1))
     if (.not.gr_dirIsAngular(JAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(JAXIS))
     end if
  end if

  if (NDIM == 3) then
     gr_minCellSizes(KAXIS) = (gr_kmax - gr_kmin) / &
          (gr_nblockZ*NZB*2**(gr_maxRefine-1))
     if (.not. gr_dirIsAngular(KAXIS)) then
        gr_minCellSize = min(gr_minCellSize,gr_minCellSizes(KAXIS))
     end if
  end if

#ifdef FL_NON_PERMANENT_GUARDCELLS
  gr_blkPtrRefCount = 0 
  gr_blkPtrRefCount_fc = 0
  gr_lastBlkPtrGotten = 0 
  gr_lastBlkPtrGotten_fc = 0
#endif

  call gr_setDataStructInfo()
  call gr_bcInit()

  !Initialize grid arrays used by IO
  allocate(gr_nToLeft(0:gr_meshNumProcs-1))
  allocate(gr_gid(nfaces+nchild+1, MAXBLOCKS))
#ifdef FLASH_GRID_PARAMESH3OR4
  allocate(gr_gsurr_blks(2,1+(K1D*2),1+(K2D*2),1+(K3D*2),MAXBLOCKS))
#endif
#ifndef BSS_GRID_ARRAYS
# if NSCRATCH_GRID_VARS > 0
  allocate(scratch(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc+1,&
       gr_kLoGc:gr_kHiGc+1,MAXBLOCKS))
# else
  allocate(scratch(1,1,1,1,1))
# endif

# if NSCRATCH_CENTER_VARS > 0
  allocate(scratch_ctr(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
# else
  allocate(scratch_ctr(1,1,1,1,1))
# endif

# if(NSCRATCH_FACEX_VARS>0)  
  allocate(scratch_facevarx( SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
# else
  allocate(scratch_facevarx(1,1,1,1,1))
# endif

# if(NSCRATCH_FACEY_VARS>0)  
  allocate(scratch_facevary( SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
       gr_kLoGc:gr_kHiGc,MAXBLOCKS))
# else
  allocate(scratch_facevary(1,1,1,1,1))
# endif  

# if(NSCRATCH_FACEZ_VARS>0)
  allocate(scratch_facevarz( SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc+K3D,MAXBLOCKS) )
# else
  allocate(scratch_facevarz(1,1,1,1,1))
# endif

  allocate(gr_xflx(NFLUXES,2,NYB,NZB,MAXBLOCKS))
  allocate(gr_yflx(NFLUXES,NXB,2,NZB,MAXBLOCKS))
  allocate(gr_zflx(NFLUXES,NXB,NYB,2,MAXBLOCKS))
  
# ifdef FLASH_HYDRO_UNSPLIT
#  if NDIM >= 2
  allocate(gr_xflx_yface(NFLUXES,2:NXB, 2   ,NZB  ,MAXBLOCKS))
  allocate(gr_yflx_xface(NFLUXES,2    ,2:NYB,NZB  ,MAXBLOCKS))
#   if NDIM == 3
  allocate(gr_xflx_zface(NFLUXES,2:NXB,NYB  , 2   ,MAXBLOCKS))
  allocate(gr_yflx_zface(NFLUXES,NXB,  2:NYB, 2   ,MAXBLOCKS))
  allocate(gr_zflx_xface(NFLUXES, 2 ,NYB    ,2:NZB,MAXBLOCKS))
  allocate(gr_zflx_yface(NFLUXES,NXB, 2     ,2:NZB,MAXBLOCKS))
#   endif
#  endif
# endif

#endif

  if(gr_meshMe == MASTER_PE) call printRefinementInfo()

contains

  subroutine printRefinementInfo()
    implicit none
    integer :: l,n
    real    :: del(MDIM)
    character(len=20) :: fmtStr
    character(len=2)  :: colHdr(MDIM) = (/'dx', 'dy', 'dz'/)

    write(*,*) 'Grid_init: resolution based on runtime params:'
    write(*,'(A9,3(A12:4x))')  'lrefine', (colHdr(n),n=1,NDIM)
    do l = lrefine_min, lrefine_max
       del (IAXIS)               = (gr_imax - gr_imin) / (gr_nblockX*NXB*2.**(l-1))
       if (NDIM > 1)  del(JAXIS) = (gr_jmax - gr_jmin) / (gr_nblockY*NYB*2.**(l-1))
       if (NDIM == 3) del(KAXIS) = (gr_kmax - gr_kmin) / (gr_nblockZ*NZB*2.**(l-1))

       if (maxval(del(IAXIS:NDIM)) .GT. 999999999999.999) then
          fmtStr = '(I7,2x,1P,3G16.3)'
       else if (minval(del(IAXIS:NDIM)) .LE. 0.0009) then
          fmtStr = '(I7,2x,1P,3G16.3)'
       else
          fmtStr = '(I7,2x,3F16.3)'
       end if

       write(*,fmtStr) l, (del(n),n=1,NDIM)
    end do
  end subroutine printRefinementInfo


end subroutine gr_initSpecific

