!!****if* source/Grid/GridMain/AMR/Amrex/gr_initNewLevelCallback
!!
!! NAME
!!  gr_initNewLevelCallback
!!
!! SYNOPSIS
!!  gr_initNewLevelCallback(integer(IN)    :: lev,
!!                          amrex_real(IN) :: time,
!!                          c_ptr(IN)      :: pba,
!!                          c_ptr(IN)      :: pdm)
!!
!! DESCRIPTION
!!  This routine is a callback routine that is registered with the AMReX AMR
!!  core at initialization.  AMReX calls this routine to create a new refinement
!!  level from scratch.
!!
!!  Specifically, for the given refinement level it
!!   (1) creates a multifab for each data type,
!!   (2) initializes these structures with the initial conditions via
!!       Simulation_initBlock,
!!   (3) fills all cell-centered guardcells, and
!!   (4) runs EoS on the interiors and GCs of all blocks to make the initial
!!       data thermodynamically consistent.
!!
!!  Note that all EoS runs are done in the mode specified by the eosModeInit
!!  runtime parameter.
!!
!!  This routine should only be invoked by AMReX.
!!
!! ARGUMENTS
!!  lev - a 0-based number identifying the refinement level to create.  The
!!        zeroth level is the coarsest level to be used in the simulation and a
!!        larger integer indicates a finer refinement.
!!  time - IGNORED
!!  pba - a C pointer to the AMReX box array object to use for constructing the
!!        multifab for the given level
!!  pdm - a C pointer to the AMReX distribution mapping of boxes across
!!        processors to be used for constructing the multifab for the given
!!        level.
!!
!! SEE ALSO
!!  Grid_initDomain
!!
!!***

#include "constants.h"
#include "Flash.h"

#if defined(FLASH_HYDRO_UNSPLIT) && defined(FLASH_UHD_HYDRO)
#include "UHD.h"
#endif

subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amr_module,          ONLY : amrex_geom
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab_build
    use amrex_fluxregister_module, ONLY : amrex_fluxregister_build
    use amrex_fillpatch_module,    ONLY : amrex_fillpatch

    use gr_physicalMultifabs,      ONLY : unk, &
                                          gr_scratchCtr, &
                                          facevarx, facevary, facevarz, &
                                          fluxes, &
                                          flux_registers
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC, &
                                          gr_preinterpolationWork, &
                                          gr_postinterpolationWork
    use flash_iterator,            ONLY : flash_iterator_t
    use flash_tile,                ONLY : flash_tile_t
    use Simulation_interface,      ONLY : Simulation_initBlock
    use Grid_interface,            ONLY : Grid_getTileIterator, &
                                          Grid_releaseTileIterator
    use Grid_data,                 ONLY : gr_eosModeInit, &
                                          gr_doFluxCorrection, &
                                          gr_interpolator, &
                                          lo_bc_amrex, hi_bc_amrex
    use Eos_interface,             ONLY : Eos_wrapped
    use Logfile_interface,         ONLY : Logfile_stamp

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm

    type(flash_iterator_t) :: itor
    type(flash_tile_t)     :: tileDesc
    real(wp), contiguous, pointer :: initData(:,:,:,:)

    integer :: n_blocks

    integer :: dir
    logical :: nodal(1:MDIM)

    integer :: i, j, k, var

    nullify(initData)

    ba = pba
    dm = pdm

    call gr_clearLevelCallback(lev)

    !!!!! CREATE MULTIFABS FOR STORING PHYSICAL DATA AT GIVEN LEVEL
    ! Cell-centered unknowns
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)

#if NFACE_VARS > 0
    ! Face variables
    nodal(:)     = .FALSE.
    nodal(IAXIS) = .TRUE.
    call amrex_multifab_build(facevarx(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#if NDIM >= 2
    nodal(:)     = .FALSE.
    nodal(JAXIS) = .TRUE.
    call amrex_multifab_build(facevary(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#if NDIM == 3
    nodal(:)     = .FALSE.
    nodal(KAXIS) = .TRUE.
    call amrex_multifab_build(facevarz(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#endif

    ! Create FABs for needed by Hydro.
    !! DEV : Control of gr_scratchCtr allocation is very hacky...
#ifdef HY_VAR2_SCRATCHCTR_VAR
# ifdef HY_XN06_SCRATCHCTR_VAR
    call amrex_multifab_build(gr_scratchCtr(lev), ba, dm, HY_XN06_SCRATCHCTR_VAR, 0)
# else
    call amrex_multifab_build(gr_scratchCtr(lev), ba, dm, HY_VAR2_SCRATCHCTR_VAR, 0)
# endif
#endif

#if NFLUXES > 0
    ! No need to store fluxes for guardcells
    do dir = 1, SIZE(fluxes, 2)
        nodal(:)   = .FALSE.
        nodal(dir) = .TRUE.
        call amrex_multifab_build(fluxes(lev, dir), ba, dm, NFLUXES, 0, &
                                  nodal=nodal)
    end do

    if ((lev > 0) .AND. (gr_doFluxCorrection)) then
        call amrex_fluxregister_build(flux_registers(lev), ba, dm, &
                                      amrex_ref_ratio(lev-1), &
                                      lev, NFLUXES)
    end if
#endif

    ! Write initial data across domain at given level and convert to conserved
    ! form where needed for proper interpolation with GC fill below.
    n_blocks = 0
    call Grid_getTileIterator(itor, ALL_BLKS, level=lev+1, tiling=.TRUE.)
    do while (itor%isValid())
        call itor%currentTile(tileDesc)

        !  We need to zero data in case we reuse blocks from previous levels
        !  but don't initialize all data in Simulation_initBlock... in particular
        !  the total vs. internal energies can cause problems in the eos call that 
        !  follows.  This includes zeroing the data in the grown tile's
        !  guardcells.
        associate(lo => tileDesc%limitsGC(LOW,  :), &
                  hi => tileDesc%limitsGC(HIGH, :))
            call tileDesc%getDataPtr(initData, CENTER)
            do           var = UNK_VARS_BEGIN, UNK_VARS_END
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            initData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do
            call tileDesc%releaseDataPtr(initData, CENTER)

#if NFACE_VARS > 0
            call tileDesc%getDataPtr(initData, FACEX)
            do           var = 1, NFACE_VARS
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)+1
                            initData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do
            call tileDesc%releaseDataPtr(initData, FACEX)
#if NDIM >= 2
            call tileDesc%getDataPtr(initData, FACEY)
            do           var = 1, NFACE_VARS
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)+1
                        do i = lo(IAXIS), hi(IAXIS)
                            initData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do
            call tileDesc%releaseDataPtr(initData, FACEY)
#endif
#if NDIM == 3
            call tileDesc%getDataPtr(initData, FACEZ)
            do           var = 1, NFACE_VARS
                do         k = lo(KAXIS), hi(KAXIS)+1
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            initData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do
            call tileDesc%releaseDataPtr(initData, FACEZ)
#endif
#endif
        end associate

        ! Give simulation the cell-centered data.  If they need to initialize
        ! face-centered data, they access it explicitly from tileDesc
        call tileDesc%getDataPtr(initData, CENTER)
        call Simulation_initBlock(initData, tileDesc)
        call tileDesc%releaseDataPtr(initData, CENTER)

        n_blocks = n_blocks + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    call Logfile_stamp(lev+1, &
          '[gr_initNewLevelCallback] Initialized data on level')

    ! Subsequent AMReX calls to gr_markRefineDerefineCallback require that the
    ! GC be filled.  We do *not* ask client code to do this, so fill GC here
    !
    ! The routine gr_estimateBlkError is only using cell-centered data for
    ! gauging refinement of blocks.   Therefore, we need not do a GC fill
    ! for the facevar[xyz] data.
    if (lev == 0) then
       ! Move all unk data to given ba/dm layout.  Do *not* use sub-cycling.
       ! -1 because of Fortran variable index starts with 1
       call amrex_fillpatch(unk(lev), time+1.0, unk(lev), &
                                      time,     unk(lev), &
                                      amrex_geom(lev), gr_fillPhysicalBC, &
                                      time, &
                                      UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)
    else
       call amrex_fillpatch(unk(lev), time+1.0, unk(lev-1), &
                                      time,     unk(lev-1), &
                                      amrex_geom(lev-1), gr_fillPhysicalBC, &
                                      time+1.0, unk(lev  ), &
                                      time,     unk(lev  ), &
                                      amrex_geom(lev  ), gr_fillPhysicalBC, &
                                      time, &
                                      UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                      amrex_ref_ratio(lev-1), &
                                      gr_interpolator, &
                                      lo_bc_amrex, hi_bc_amrex, &
                                      gr_preinterpolationWork, &
                                      gr_postinterpolationWork)
    end if

    call Logfile_stamp(lev+1, &
          '[gr_initNewLevelCallback] GC fill')

    ! Run EoS on interiors and GCs in preparation for refinement check
    call Grid_getTileIterator(itor, ALL_BLKS, level=lev+1, tiling=.TRUE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       call tileDesc%getDataPtr(initData, CENTER)
       call Eos_wrapped(gr_eosModeInit, tileDesc%limitsGC, initData)
       call tileDesc%releaseDataPtr(initData, CENTER)

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    write(*,'(A,I10,A,I0)') "Created and initialized ", n_blocks, &
                           " blocks on level ", lev + 1

end subroutine gr_initNewLevelCallback

