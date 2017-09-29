#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"
#include "constants.h"

subroutine gr_markRefineDerefineCallback(lev, tags, time, tagval, clearval) bind(c)
   use iso_c_binding
   use amrex_fort_module,      ONLY : wp => amrex_real
   use amrex_box_module,       ONLY : amrex_box
   use amrex_tagbox_module,    ONLY : amrex_tagboxarray
   use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy

   use Grid_data,              ONLY : gr_numRefineVars, &
                                      gr_refine_cutoff, gr_derefine_cutoff, &
                                      gr_refine_filter, gr_refine_var, &
                                      gr_maxRefine, gr_enforceMaxRefinement
   use Grid_interface,         ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
   use gr_interface,           ONLY : gr_estimateBlkError
   use gr_physicalMultifabs,   ONLY : unk
   use block_metadata,         ONLY : block_metadata_t

   implicit none
 
   integer,           intent(IN), value :: lev
   type(c_ptr),       intent(IN), value :: tags 
   real(wp),          intent(IN), value :: time
   character(c_char), intent(IN), value :: tagval
   character(c_char), intent(IN), value :: clearval

   type(amrex_tagboxarray) :: tag
   type(amrex_mfiter)      :: mfi                                                             
   type(amrex_box)         :: bx
   type(block_metadata_t)  :: blockDesc

   real(wp),               contiguous, pointer :: solnData(:,:,:,:)
   character(kind=c_char), contiguous, pointer :: tagData(:,:,:,:)

   real, allocatable :: errors(:)
   real              :: refineCut, derefineCut, refineFilter

   integer :: off(1:MDIM)

   integer :: iref
   integer :: i, j, k, l

#ifdef DEBUG_GRID
   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Started on level ", lev + 1
#endif
 
   tag = tags

   allocate(errors(gr_numRefineVars))

   !DEVNOTE:  Can test with tiling later - KW
   ! unk used the same 0-based level indexing used here by AMReX
   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)
   do while(mfi%next())
      bx = mfi%tilebox()

      ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
      ! Level must be 1-based index and limits/limitsGC must be 1-based also
      ! DEVNOTE: Should we use gr_[ijk]guard here?
      blockDesc%level = lev + 1
      blockDesc%grid_index = mfi%grid_index()
      blockDesc%limits(LOW,  :) = 1
      blockDesc%limits(HIGH, :) = 1
      blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
      blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
      blockDesc%limitsGC(LOW,  :) = 1
      blockDesc%limitsGC(HIGH, :) = 1
      blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
      blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

      errors(:) = 0.0d0
      do l = 1, gr_numRefineVars
         iref = gr_refine_var(l)
         refineFilter = gr_refine_filter(l)
         call gr_estimateBlkError(errors(l), blockDesc, iref, refineFilter)
      end do

      call Grid_getBlkPtr(blockDesc, solnData, CENTER)

      associate (lo   => blockDesc%limits(LOW,  :), &
                 hi   => blockDesc%limits(HIGH, :), &
                 loGC => blockDesc%limitsGC(LOW,  :), &
                 hiGC => blockDesc%limitsGC(HIGH, :))

        ! tagData is one cell larger on all borders than interior and 0-based
        ! Shift to 1-based here
        off = lo
        off(1:NDIM) = lo(1:NDIM) - 1
        tagData(off(1):, off(2):, off(3):, 1:) => tag%dataptr(mfi)

#ifdef DEBUG_TAGDATA
        print*,'markRD_cb: lbound(solnData):', lbound(solnData)
        print*,'markRD_cb: ubound(solnData):', ubound(solnData)
        print*,'markRD_cb: lbound(tagData):', lbound(tagData)
        print*,'markRD_cb: ubound(tagData):', ubound(tagData)
        print*,'markRD_cb: tagData in  =', tagData
#endif

        tagData(:, :, :, :) = clearval
 rloop: do l = 1, gr_numRefineVars
            if (gr_refine_var(l) < 1)   CYCLE

            ! DEV: Not clear if we can tag cells to inform AMReX that
            ! we request derefinement.  TODO Figure out.
            if (errors(l) > gr_refine_cutoff(l)) then
                ! Tag single cell in block that is not on boundary
                i = INT(0.5d0 * DBLE(lo(IAXIS) + hi(IAXIS)))
                j = INT(0.5d0 * DBLE(lo(JAXIS) + hi(JAXIS)))
                k = INT(0.5d0 * DBLE(lo(KAXIS) + hi(KAXIS)))

                ! NOTE: last dimension has range 1:1
                tagData(i, j, k, 1) = tagval

#ifdef DEBUG_GRID
                write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                                    "      Tag block for refinement at level", &
                                    (lev+1)
                write(*,'(A,A,I4,I4,I4)') "[gr_markRefineDerefineCallback]", &
                                          "         lower: ", &
                                          lo(IAXIS), lo(JAXIS), lo(KAXIS)
                write(*,'(A,A,I4,I4,I4)') "[gr_markRefineDerefineCallback]", &
                                          "         upper: ", &
                                          hi(IAXIS), hi(JAXIS), hi(KAXIS)
                write(*,'(A,A,I4,I4,I4)') "[gr_markRefineDerefineCallback]", &
                                          "         Tag cell ", i, j, k 
#endif
                
                EXIT rloop
            end if
        end do rloop
      end associate

      call Grid_releaseBlkPtr(blockDesc, solnData)
   end do
   call amrex_mfiter_destroy(mfi)

   deallocate(errors)

#ifdef DEBUG_GRID
   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
#endif

end subroutine gr_markRefineDerefineCallback

