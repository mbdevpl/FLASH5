!!****if* source/Grid/GridMain/AMR/Amrex/gr_markRefineDerefineCallback
!!
!! NAME
!!  gr_markRefineDerefineCallback
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefineCallback(integer(IN) :: lev,
!!                                c_ptr(IN)   :: tags,
!!                                real(IN)    :: time,
!!                                c_char(IN)  :: tagval,
!!                                c_char(IN)  :: clearval)
!!  
!!  DESCRIPTION
!!  
!!  This routine is a callback subroutine that is registered with AMReX's
!!  AMR Core layer at initialization.  AMReX may call this subroutine many times
!!  during the process of grid refinement so that FLASH may communicate which
!!  blocks in the given level require refinement.  The final refinement
!   decisions are made by AMReX based on the information gathered with this
!!  callback.
!!
!!  This routine iterates across all blocks in the given level and determines if
!!  the current block needs refinement.  If it does, then all cells in the AMReX
!!  tagbox associated with the block are marked for refinement by setting their
!!  value to tagval.  If not, then all cells are set to clearval.
!!
!!  A block is marked for refinement if the block's error estimate for any
!!  refinement variable is greater than the variable's associated refinement
!!  cutoff value.
!!
!!  ARGUMENTS 
!!
!!    lev - the 0-based level index
!!    tags - C-pointer to an AMReX tagbox array.  The elements of this are tag
!!           boxes.  The cells of these tagboxes are set to communicate a need
!!           to refine the associated block.
!!    time - not used with FLASH
!!    tagval - for full, rich AMReX tagging, this values should be assigned to
!!             each cell that has insufficient resolution.
!!    clearval - for full, rich AMReX tagging, this values should be assigned to
!!               each cell that has sufficient resolution.
!! 
!!  SEE ALSO
!!  
!!    gr_estimateBlkError
!!
!!***

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
                                      gr_maxRefine, gr_enforceMaxRefinement, &
                                      gr_minRefine
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

   real :: error
   real :: refineCut, derefineCut, refineFilter

   integer :: off(1:MDIM)

   integer :: iref
   integer :: i, j, k, l

#ifdef DEBUG_GRID
   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Started on level ", lev + 1
#endif
 
   tag = tags

   if (lev < gr_minRefine - 1) then
#ifdef DEBUG_GRID
      write(*,'(A,A,I4,I4,I4)') "[gr_markRefineDerefineCallback]", &
                                "         derefinement to this level not allowed"
#endif

      ! Enforce 1-based minimum level contraint
      call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)

      do while(mfi%next())
         tagData => tag%dataptr(mfi)
         tagData(:, :, :, :) = tagval
         nullify(tagData)
      end do

      call amrex_mfiter_destroy(mfi)
      RETURN
   end if

   !DEVNOTE:  Can test with tiling later - KW
   ! unk used the same 0-based level indexing used here by AMReX
   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)
   do while(mfi%next())
      bx = mfi%fabbox()

      ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
      ! Level must be 1-based index and limits/limitsGC must be 1-based also
      blockDesc%level = lev + 1
      blockDesc%grid_index = mfi%grid_index()
      blockDesc%limits(LOW,  :) = 1
      blockDesc%limits(HIGH, :) = 1
      blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1 + NGUARD
      blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1 - NGUARD
      blockDesc%limitsGC(LOW,  :) = 1
      blockDesc%limitsGC(HIGH, :) = 1
      blockDesc%limitsGC(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
      blockDesc%limitsGC(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1

      call Grid_getBlkPtr(blockDesc, solnData, CENTER)

      ! The width of the halo of gaurdcells included in the tagbox array may
      ! differ.  According to Weiqun, this is necessary for ensuring proper
      ! nesting.
      tagData => tag%dataptr(mfi)
      
      associate (lo     => blockDesc%limits(LOW,  :), &
                 hi     => blockDesc%limits(HIGH, :))

#ifdef DEBUG_TAGDATA
        print*,'markRD_cb: lbound(solnData):', lbound(solnData)
        print*,'markRD_cb: ubound(solnData):', ubound(solnData)
        print*,'markRD_cb: lbound(tagData):', (lbound(tagData) + 1)
        print*,'markRD_cb: ubound(tagData):', (ubound(tagData) + 1)
#endif

        tagData(:, :, :, :) = clearval
 rloop: do l = 1, gr_numRefineVars
            iref = gr_refine_var(l)
            if (iref < 1)   CYCLE
    
            error = 0.0d0
            refineFilter = gr_refine_filter(l)
            call gr_estimateBlkError(error, blockDesc, iref, refineFilter)

            ! Refinement is based on Berger-Rigoutsis algorithm, for which each
            ! cell is marked as having sufficient or insufficient resolution.
            ! There is no means to indicate derefine/stay/refine as with
            ! Paramesh.
            if (error > gr_refine_cutoff(l)) then
                ! According to Weiqun, when AMReX is setup in octree mode,
                ! tagging a single cell in a block is sufficient for indicating
                ! a need to refine.
                !
                ! We err on the side of caution by tagging all cells in box
                ! array to ensure octree refinement
                tagData(:, :, :, :) = tagval

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

        nullify(tagData)
      end associate

      call Grid_releaseBlkPtr(blockDesc, solnData)
   end do
   call amrex_mfiter_destroy(mfi)

#ifdef DEBUG_GRID
   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
#endif

end subroutine gr_markRefineDerefineCallback

