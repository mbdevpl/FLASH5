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

   implicit none
 
#include "Flash.h"
#include "constants.h"

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

   integer :: off(1:MDIM)

   integer :: refine_to
   integer :: i, j

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Started on level ", lev + 1
   
   tag = tags

   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)                             
   do while(mfi%next())
      bx = mfi%tilebox()

      ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
      ! Level must be 1-based index and limits/limitsGC must be 1-based also
      ! DEVNOTE: Should we use gr_[ijk]guard here?
      blockDesc%level = lev + 1
      blockDesc%grid_index = -1
      blockDesc%limits(LOW,  :) = 1
      blockDesc%limits(HIGH, :) = 1
      blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
      blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
      blockDesc%limitsGC(LOW,  :) = 1
      blockDesc%limitsGC(HIGH, :) = 1
      blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
      blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD
 
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

        ! Force refinement so that one corner is refined more than others
        tagData(:, :, :, :) = clearval
        if (      (lev < 2) &
            .AND. (lo(IAXIS) <= 2) .AND. (2 <= hi(IAXIS)) &
            .AND. (lo(JAXIS) <= 2) .AND. (2 <= hi(JAXIS))) then
            tagData(2, 2, 1, 1) = tagval
        end if
      end associate
   end do
   call amrex_mfiter_destroy(mfi)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

