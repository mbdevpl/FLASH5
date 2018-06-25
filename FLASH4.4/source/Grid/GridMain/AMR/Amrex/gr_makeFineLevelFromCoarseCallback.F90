!!****if* source/Grid/GridMain/AMR/Amrex/gr_makeFineLevelFromCoarseCallback
!!
!! NAME
!!
!!  gr_makeFineLevelFromCoarseCallback
!!
!! SYNOPSIS
!!
!!  gr_makeFineLevelFromCoarseCallback(integer(IN)    :: lev,
!!                                     amrex_real(IN) :: time,
!!                                     c_ptr(IN)      :: pba,
!!                                     c_ptr(IN)      :: pdm)
!!
!! DESCRIPTION
!!
!!  This routine is a callback routine that is registered with the AMReX AMR
!!  core at initialization.  AMReX calls this routine to populate a
!!  newly-created level with data prolongated from parent blocks and possibly
!!  with boundary data.  Prolongation is accomplished with AMReX's conservative
!!  linear interpolation routine.
!!
!!  It is assumed that, where applicable, the data is in conserved form.
!!  Upon returning, the remade multifab will have data in all interiors as
!!  well as guardcells.  Note that since the data is still in conserved form
!!  and EoS might expect primitive form, EoS is not run.  It is therefore the
!!  responsibility of the caller to manage this step.
!!
!!  NOTE: This implementation, while presently functional, is incorrect.  A user
!!        could supply their own BC routine that sensibly assumes that data is
!!        in primitive form.  However, here we pass the BC routines data in
!!        conserved form.
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

#include "constants.h"
#include "Flash.h"

#if defined(FLASH_HYDRO_UNSPLIT) && defined(FLASH_UHD_HYDRO)
#include "UHD.h"
#endif

subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab_build, &
                                          amrex_mfiter, &
                                          amrex_mfiter_build, &
                                          amrex_mfiter_destroy
    use amrex_fluxregister_module, ONLY : amrex_fluxregister_build
    use amrex_fillpatch_module,    ONLY : amrex_fillcoarsepatch
    use amrex_interpolater_module, ONLY : amrex_interp_cell_cons

    use Grid_data,                 ONLY : gr_doFluxCorrection, &
                                          gr_amrexDidRefinement, &
                                          lo_bc_amrex, hi_bc_amrex
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk, &
                                          gr_scratchCtr, &
                                          facevarx, facevary, facevarz, &
                                          fluxes, &
                                          flux_registers

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter)    :: mfi

    integer :: nFab

    integer :: dir
    logical :: nodal(1:MDIM)

    ! Communicate to Grid_updateRefinement that we have created a level
    gr_amrexDidRefinement = .TRUE.

    ba = pba
    dm = pdm

    !!!!!----- (Re)create FABS for storing physical data at this level
    call gr_clearLevelCallback(lev)

    call amrex_multifab_build(unk     (lev), ba, dm, NUNK_VARS, NGUARD)
#if NFACE_VARS > 0
    ! DEVNOTE: TODO Create these wrt proper face-centered boxes
    call amrex_multifab_build(facevarx(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevary(lev), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarz(lev), ba, dm, NUNK_VARS, NGUARD)
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

    !!!!!----- Fill new refinement level via interpolation from parent block
    ! This *hopefully* will do the guard cell fill as well
    ! NOTE: FLASH does not use sub-cycling (temporal interpolation)
    !
    ! -1 because Fortran variable index starts with 1
    call amrex_fillcoarsepatch(unk(lev), time,     unk(lev-1),  &
                                         time+0.1, unk(lev-1),  &
                                         amrex_geom(lev-1), gr_fillPhysicalBC,  &
                                         amrex_geom(lev  ), gr_fillPhysicalBC,  &
                                         time, &
                                         UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                         amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                         lo_bc_amrex, hi_bc_amrex) 

    nFab = 0
    call amrex_mfiter_build(mfi, unk(lev), tiling=.false.)
    do while(mfi%next())
        nFab = nFab + 1 
    end do
    call amrex_mfiter_destroy(mfi)

    write(*,'(A,I0,A,I0,A)') "Made fine level ", lev + 1, " - ", nFab, " blocks"

end subroutine gr_makeFineLevelFromCoarseCallback

