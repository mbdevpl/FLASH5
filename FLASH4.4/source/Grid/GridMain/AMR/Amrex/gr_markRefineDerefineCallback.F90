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

   integer :: i, j, k

   tag = tags
   
   if (time == 0.0d0)  RETURN

   call amrex_mfiter_build(mfi, unk(lev), tiling=.true.)                             
   do while(mfi%next())
      bx = mfi%tilebox()
      solnData => unk(lev)%dataptr(mfi)
      tagData  => tag%dataptr(mfi)

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
                 hi => blockDesc%limits(HIGH, :))
        tagData(:, :, :, 1) = clearval
        do       k = lo(3), hi(3)
           do    j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 if (solnData(i, j, k, 1) >= 1.1_wp) then
                    tagData(i, j, k, 1) = tagval
                 end if
              end do
           end do
        end do
      end associate
   end do
   call amrex_mfiter_destroy(mfi)

   write(*,*) "[gr_markRefineDerefineCallback] on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

