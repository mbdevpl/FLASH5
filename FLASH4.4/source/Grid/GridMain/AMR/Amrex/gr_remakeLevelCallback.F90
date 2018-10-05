!!****if* source/Grid/GridMain/AMR/Amrex/gr_remakeLevelCallback
!!
!! NAME
!!  gr_remakeLevelCallback
!!
!! SYNOPSIS
!!  gr_remakeLevelCallback(integer(IN)    :: lev,
!!                         amrex_real(IN) :: time,
!!                         c_ptr(IN)      :: pba,
!!                         c_ptr(IN)      :: pdm)
!!
!! DESCRIPTION
!!  This routine is a callback routine that is registered with the AMReX AMR
!!  core at initialization.  AMReX calls this routine to reestablish the data in
!!  a multifab at the given level onto a new multifab specified through the given
!!  box array and distribution map.
!!
!!  Upon returning, the remade multifab will have data in all interiors 
!!  as well as guardcells.  EoS is not run on blocks so that it is the 
!!  responsibility of the code that triggered this callback to execute
!!  EoS where needed.
!!
!!  In detail, for the given refinement level this routine
!!   (1) uses AMReX patchfill to copy data from each original multifab to an
!!       associated buffer multifab built with the new box layout and 
!!       distribution mapping,
!!   (2) rebuild each original multifab, and
!!   (3) copy the new interior/GC data from the buffer to the associated 
!!       rebuilt multifab.
!!
!!  Note that step (1) might require that AMReX execute prolongation operations
!!  using the AMReX conservative linear interpolation algorithm if new boxes are
!!  added to the level.  
!!
!!  The AMReX pathcfill routines called by this routine might trigger calls
!!  to gr_fillPhysicalBC, which requires that data be in primitive form.
!!  Therefore, this routine assumes that no primitive-to-conservative
!!  transformations have been done.  AMReX will use a FLASH-provided routine to
!!  convert primitive form data to conservative form data before interpolation.
!!  It will use a different FLASH-provided routine to convert interpolated
!!  conservative form data to primitive form data as well.
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
!!  gr_preinterpolationWork
!!  gr_postinterpolationWork
!!
!!***

#include "Flash.h"
#include "constants.h"

#if defined(FLASH_HYDRO_UNSPLIT) && defined(FLASH_UHD_HYDRO)
#include "UHD.h"
#endif

subroutine gr_remakeLevelCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_box_module,          ONLY : amrex_box
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab, &
                                          amrex_multifab_build, &
                                          amrex_multifab_destroy, &
                                          amrex_mfiter, &
                                          amrex_mfiter_build, &
                                          amrex_mfiter_destroy
    use amrex_fluxregister_module, ONLY : amrex_fluxregister_build
    use amrex_fillpatch_module,    ONLY : amrex_fillpatch

    use Grid_data,                 ONLY : lo_bc_amrex, hi_bc_amrex, &
                                          gr_amrexDidRefinement, &
                                          gr_doFluxCorrection, &
                                          gr_interpolator
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_preinterpolationWork, &
                                          gr_postinterpolationWork, &
                                          gr_fillPhysicalFaceBC, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz, &
                                          gr_scratchCtr, &
                                          fluxes, &
                                          flux_registers

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_box)       :: bx
    type(amrex_multifab)  :: tmp_unk
#if NFACE_VARS > 0
    type(amrex_multifab)  :: tmp_facevarx
#if NDIM >= 2
    type(amrex_multifab)  :: tmp_facevary
#endif
#if NDIM == 3
    type(amrex_multifab)  :: tmp_facevarz
#endif
#endif

    type(amrex_mfiter) :: mfi
    integer            :: nFab

    integer :: dir
    logical :: nodal(1:MDIM)

    ! Communicate to Grid_updateRefinement that we are regridding
    gr_amrexDidRefinement = .TRUE.

    ba = pba
    dm = pdm

    !!!!! SAVE DATA IN BUFFER WITH GIVEN BOXARRAY/DISTRIBUTION
    ! Get all unk interior data
    call amrex_multifab_build(tmp_unk, ba, dm, NUNK_VARS, NGUARD)

#if NFACE_VARS > 0
    ! Face variables
    nodal(:)     = .FALSE.
    nodal(IAXIS) = .TRUE.
    call amrex_multifab_build(tmp_facevarx, ba, dm, NFACE_VARS, NGUARD, nodal)
#if NDIM >= 2
    nodal(:)     = .FALSE.
    nodal(JAXIS) = .TRUE.
    call amrex_multifab_build(tmp_facevary, ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#if NDIM == 3
    nodal(:)     = .FALSE.
    nodal(KAXIS) = .TRUE.
    call amrex_multifab_build(tmp_facevarz, ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#endif

    if (lev == 0) then
       ! Move all unk data (interior and GC) to given ba/dm layout.
       ! Do *not* use sub-cycling.
       call amrex_fillpatch(tmp_unk, time+1.0d0, unk(lev), &
                                     time,       unk(lev), &
                                     amrex_geom(lev), gr_fillPhysicalBC, &
                                     time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)       
#if NFACE_VARS > 0
       call amrex_fillpatch(tmp_facevarx, time+1.0, facevarx(lev), &
                                          time,     facevarx(lev), &
                                          amrex_geom(lev), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS)       
#if NDIM >= 2
       call amrex_fillpatch(tmp_facevary, time+1.0, facevary(lev), &
                                          time,     facevary(lev), &
                                          amrex_geom(lev), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS)       
#endif
#if NDIM == 3
       call amrex_fillpatch(tmp_facevarz, time+1.0, facevarz(lev), &
                                          time,     facevarz(lev), &
                                          amrex_geom(lev), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS)       
#endif
#endif
    else
       call amrex_fillpatch(tmp_unk, time+1.0, unk(lev-1), &
                                     time,     unk(lev-1), &
                                     amrex_geom(lev-1), gr_fillPhysicalBC, &
                                     time+1.0, unk(lev  ), &
                                     time,     unk(lev  ), &
                                     amrex_geom(lev  ), gr_fillPhysicalBC, &
                                     time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                     amrex_ref_ratio(lev-1), gr_interpolator, &
                                     lo_bc_amrex, hi_bc_amrex, &
                                     gr_preinterpolationWork, &
                                     gr_postinterpolationWork)

#if NFACE_VARS > 0
       call amrex_fillpatch(tmp_facevarx, time+1.0, facevarx(lev-1), &
                                          time,     facevarx(lev-1), &
                                          amrex_geom(lev-1), gr_fillPhysicalBC, &
                                          time+1.0, facevarx(lev  ), &
                                          time,     facevarx(lev  ), &
                                          amrex_geom(lev  ), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS, &
                                          amrex_ref_ratio(lev-1), gr_interpolator, &
                                          lo_bc_amrex, hi_bc_amrex)       
#if NDIM >= 2
       call amrex_fillpatch(tmp_facevary, time+1.0, facevary(lev-1), &
                                          time,     facevary(lev-1), &
                                          amrex_geom(lev-1), gr_fillPhysicalBC, &
                                          time+1.0, facevary(lev  ), &
                                          time,     facevary(lev  ), &
                                          amrex_geom(lev  ), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS, &
                                          amrex_ref_ratio(lev-1), gr_interpolator, &
                                          lo_bc_amrex, hi_bc_amrex)       
#endif
#if NDIM == 3
       call amrex_fillpatch(tmp_facevarz, time+1.0, facevarz(lev-1), &
                                          time,     facevarz(lev-1), &
                                          amrex_geom(lev-1), gr_fillPhysicalBC, &
                                          time+1.0, facevarz(lev  ), &
                                          time,     facevarz(lev  ), &
                                          amrex_geom(lev  ), gr_fillPhysicalBC, &
                                          time, 1, 1, NFACE_VARS, &
                                          amrex_ref_ratio(lev-1), gr_interpolator, &
                                          lo_bc_amrex, hi_bc_amrex)       
#endif
#endif
    end if

    !!!!! REBUILD MFABS AT LEVEL AND FILL FROM BUFFERS
    ! If GC are not copied, then Hydro_computeDt fails
    call gr_clearLevelCallback(lev)

    ! cell-centered data
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)
    call unk(lev)%copy(tmp_unk, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, NGUARD)
    call amrex_multifab_destroy(tmp_unk)

#if NFACE_VARS > 0
    ! Face-centered data
    nodal(:)     = .FALSE.
    nodal(IAXIS) = .TRUE.
    call amrex_multifab_build(facevarx(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
    call facevarx(lev)%copy(tmp_facevarx, 1, 1, NFACE_VARS, NGUARD)
    call amrex_multifab_destroy(tmp_facevarx)
#if NDIM >= 2
    nodal(:)     = .FALSE.
    nodal(JAXIS) = .TRUE.
    call amrex_multifab_build(facevary(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
    call facevary(lev)%copy(tmp_facevary, 1, 1, NFACE_VARS, NGUARD)
    call amrex_multifab_destroy(tmp_facevary)
#endif
#if NDIM == 3
    nodal(:)     = .FALSE.
    nodal(KAXIS) = .TRUE.
    call amrex_multifab_build(facevarz(lev), ba, dm, NFACE_VARS, NGUARD, nodal)
    call facevarz(lev)%copy(tmp_facevarz, 1, 1, NFACE_VARS, NGUARD)
    call amrex_multifab_destroy(tmp_facevarz)
#endif
#endif

    nFab = 0
    call amrex_mfiter_build(mfi, unk(lev), tiling=.false.)
    do while(mfi%next())
        nFab = nFab + 1 
    end do
    call amrex_mfiter_destroy(mfi)

    !! DEV : Control of gr_scratchCtr allocation is very hacky...
#ifdef HY_VAR2_SCRATCHCTR_VAR
# ifdef HY_XN06_SCRATCHCTR_VAR
    call amrex_multifab_build(gr_scratchCtr(lev), ba, dm, HY_XN06_SCRATCHCTR_VAR, 0)
# else
    call amrex_multifab_build(gr_scratchCtr(lev), ba, dm, HY_VAR2_SCRATCHCTR_VAR, 0)
# endif
#endif

#if NFLUXES > 0
    !!!!! REBUILD FLUX MFABS WITHOUT UPDATING DATA
    ! No need to store fluxes for guardcells
    do dir = 1, SIZE(fluxes, 2)
        nodal(:)   = .FALSE.
        nodal(dir) = .TRUE.
        call amrex_multifab_build(fluxes(lev, dir), ba, dm, NFLUXES, 0, &
                                  nodal=nodal)
    end do

    !!!!! REBUILD FLUX REGISTER
    if ((lev > 0) .AND. (gr_doFluxCorrection)) then
        call amrex_fluxregister_build(flux_registers(lev), ba, dm, &
                                      amrex_ref_ratio(lev-1), &
                                      lev, NFLUXES)
    end if
#endif

    write(*,'(A,I0,A,I0,A)') "Remade level ", (lev+1), " - ", nFab, " blocks"

end subroutine gr_remakeLevelCallback

