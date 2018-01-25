!!****if* source/Grid/GridMain/AMR/Amrex/gr_remakeLevelCallback
!!
!! NAME
!!
!!  gr_remakeLevelCallback
!!
!! SYNOPSIS
!!
!!  gr_remakeLevelCallback(integer(IN)    :: lev,
!!                         amrex_real(IN) :: time,
!!                         c_ptr(IN)      :: pba,
!!                         c_ptr(IN)      :: pdm)
!!
!! DESCRIPTION
!!
!!  This routine is a callback routine that is registered with the AMReX AMR
!!  core at initialization.  AMReX calls this routine to reestablish the data in
!!  a multifab at the given level onto a new multifab specified through the given
!!  box array and distribution map.
!!
!!  It is assumed that the multifab data at the given level is correct and has
!!  had EoS run on it.  Upon returning, the remade multifab will have correct
!!  data on all interiors/guardcells and these will have had EoS run on them.
!!
!!  In detail, for the given refinement level this routine
!!   (1) uses AMReX patchfill to copy data from the original multifab to a
!!       buffer multifab built with the new box layout and distribution mapping,
!!   (2) rebuild the original multifab,
!!   (3) copy the new interior/GC data from the buffer to the rebuilt multifab, and
!!   (4) run EoS on the interiors of all blocks to make the data 
!!       thermodynamically consistent.
!!
!!  Note that step (1) might require that AMReX execute prolongation operations
!!  using the AMReX conservative linear interpolation algorithm if new boxes are
!!  added to the level.  All EoS runs are done in the mode specified by the
!!  eosMode runtime parameter.
!!
!!  This routine should only be invoked by AMReX.
!!
!! ARGUMENTS
!!
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
!!***

#include "Flash.h"
#include "constants.h"

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
    use amrex_fillpatch_module,    ONLY : amrex_fillpatch
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use Grid_data,                 ONLY : lo_bc_amrex, hi_bc_amrex, &
                                          gr_eosMode, &
                                          gr_amrexDidRefinement
    use Grid_interface,            ONLY : Grid_getBlkIterator, &
                                          Grid_releaseBlkIterator, &
                                          Grid_getBlkPtr, Grid_releaseBlkPtr
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk
    use block_iterator,            ONLY : block_iterator_t
    use block_metadata,            ONLY : block_metadata_t
    use Eos_interface,             ONLY : Eos_wrapped

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_box)       :: bx
    type(amrex_multifab)  :: mfab

    type(block_iterator_t)        :: itor
    type(block_metadata_t)        :: blockDesc
    real(wp), contiguous, pointer :: solnData(:,:,:,:)
    integer                       :: nFab

    ! Communicate to Grid_updateRefinement that we are regridding
    gr_amrexDidRefinement = .TRUE.

    ba = pba
    dm = pdm

    !!!!! SAVE DATA IN BUFFER WITH GIVEN BOXARRAY/DISTRIBUTION
    ! Get all unk interior data
    call amrex_multifab_build(mfab, ba, dm, NUNK_VARS, NGUARD)
    ! DEVNOTE: TODO Include facevars in this process

    if (lev == 0) then
       ! Move all unk data (interior and GC) to given ba/dm layout.
       ! Do *not* use sub-cycling.
       call amrex_fillpatch(mfab, time+1.0d0, unk(lev), &
                                  time,       unk(lev), &
                                  amrex_geom(lev), gr_fillPhysicalBC, &
                                  time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)       
    else
       call amrex_fillpatch(mfab, time+1.0d0, unk(lev-1), &
                                  time,       unk(lev-1), &
                                  amrex_geom(lev-1), gr_fillPhysicalBC, &
                                  time+1.0e0, unk(lev  ), &
                                  time,       unk(lev  ), &
                                  amrex_geom(lev  ), gr_fillPhysicalBC, &
                                  time, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                  amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                  lo_bc_amrex, hi_bc_amrex)       
    end if

    !!!!! REBUILD MFAB AT LEVEL AND FILL FROM BUFFER
    call gr_clearLevelCallback(lev)
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)

    ! If GC are not copied, then Hydro_computeDt fails
    call unk(lev)%copy(mfab, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, NGUARD)

    call amrex_multifab_destroy(mfab)

    nFab = 0
    call Grid_getBlkIterator(itor, ALL_BLKS, level=lev+1, tiling=.FALSE.)
    do while (itor%is_valid())
       call itor%blkMetaData(blockDesc)

       call Grid_getBlkPtr(blockDesc, solnData, CENTER)
       call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)
       call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)

        nFab = nFab + 1 
       call itor%next()
    end do
    call Grid_releaseBlkIterator(itor)

    write(*,'(A,I0,A,I0,A)') "Remade level ", (lev+1), " - ", nFab, " blocks"

end subroutine gr_remakeLevelCallback

