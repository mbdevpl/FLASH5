subroutine gr_markRefineDerefineCallback(lev, tags, time, tagval, clearval) bind(c)
   use iso_c_binding
   use amrex_fort_module,      ONLY : wp => amrex_real
   use amrex_box_module,       ONLY : amrex_box
   use amrex_tagbox_module,    ONLY : amrex_tagboxarray
   use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy, &
                                      amrex_multifab_build
 
   use block_metadata,         ONLY : block_metadata_t
   use gr_physicalMultifabs,   ONLY : unk
  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine
   use gr_interface,           ONLY : gr_estimateBlkError  ! to be used RSN

   implicit none
 
#include "constants.h"
#include "Flash.h"

   integer,           intent(IN), value :: lev
   type(c_ptr),       intent(in), value :: tags 
   real(wp),          intent(in), value :: time
   character(c_char), intent(in), value :: tagval
   character(c_char), intent(in), value :: clearval

   type(amrex_tagboxarray) :: tag
   type(amrex_mfiter)      :: mfi                                                             
   type(amrex_box)         :: bx
   type(block_metadata_t)  :: blockDesc

   real(wp),               contiguous, pointer :: solnData(:,:,:,:)
   character(kind=c_char), contiguous, pointer :: tagData(:,:,:,:)

   real, allocatable :: errors(:)
   real              :: refineCut, derefineCut, refineFilter

   integer :: off(1:MDIM)

   integer :: i, j, k, iref
   integer :: l
   logical :: refine, derefine, stay

   tag = tags
   
!!$   if (time == 0.0d0)  RETURN     ! DEV: Why? - KW

   allocate(errors(gr_numRefineVars))

   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.) !DEVNOTE:  Can test with tiling later - KW
   do while(mfi%next())
      bx = mfi%tilebox()

      ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
      ! Level must be 1-based index and limits/limitsGC must be 1-based also
      ! DEVNOTE: Should we use gr_[ijk]guard here?
      blockDesc%level = lev + 1
      blockDesc%grid_index = mfi%grid_index() ! DEVNOTE: We can use this in AMREX or AMREXTRANSITION mode
      blockDesc%limits(LOW,  :) = 1
      blockDesc%limits(HIGH, :) = 1
      blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
      blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
      blockDesc%limitsGC(LOW,  :) = 1
      blockDesc%limitsGC(HIGH, :) = 1
      blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
      blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

      do l = 1,gr_numRefineVars
         iref = gr_refine_var(l)
         refineFilter = gr_refine_filter(l)
         call gr_estimateBlkError(errors(l),blockDesc,iref,refineFilter)
      end do
 
      associate (lo => blockDesc%limits(LOW,  :), &
                 hi => blockDesc%limits(HIGH, :), &
                 loGC => blockDesc%limitsGC(LOW,  :), &
                 hiGC => blockDesc%limitsGC(HIGH, :))
        ! Makes this 1-based cell indexing
        solnData(loGC(1):, loGC(2):, loGC(3):, 1:) => unk(lev)%dataptr(mfi)

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
        refine   = .FALSE.
        derefine = .FALSE.
        stay     = .FALSE.
        do l = 1,gr_numRefineVars
           iref = gr_refine_var(l)
           if (iref < 1) CYCLE
           refineCut = gr_refine_cutoff(l)
           derefineCut = gr_derefine_cutoff(l)

           if (.not.refine.and. .not.stay &
                &          .and.(errors(l) .le. derefineCut)) then
              derefine = .TRUE.
           else
              derefine = .FALSE.
           end if

           if (errors(l) > refineCut) then
              derefine = .FALSE.
              refine = .TRUE.
           end if

           if (errors(l) .gt. derefineCut)  &
                &           stay = .TRUE.

           ! DEV: The following relies on gr_maxRefine being properly set, which appears to be not
           ! always the case.
#if(0)
           if (gr_enforceMaxRefinement .AND. blockDesc%level.ge.gr_maxRefine)  &
                &           refine = .FALSE.
#endif



        end do
        if (refine) then
           do         k = lo(3), hi(3)
              do     j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    tagData(i, j, k, 1) = tagval ! Note: last dimension has the range 1:1
                 end do
              end do
           end do
        else if (derefine) then
           write(*,*) "Untag (how??) at index",blockDesc%grid_index,", ref,deref,stay is", refine,derefine,stay
        end if
#ifdef DEBUG_TAGDATA
        print*,'markRD_cb: tagData out =', tagData
#endif
      end associate
   end do
   call amrex_mfiter_destroy(mfi)

   deallocate(errors)

   write(*,*) "[gr_markRefineDerefineCallback] Finished on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

