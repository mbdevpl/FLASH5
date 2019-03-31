#include "Flash.h"
#include "constants.h"

subroutine gr_markRefineDerefineCallback(lev, tags, time, tagval, clearval) bind(c)
   use iso_c_binding
   use amrex_fort_module,     ONLY : wp => amrex_real
   use amrex_box_module,      ONLY : amrex_box
   use amrex_tagbox_module,   ONLY : amrex_tagboxarray
   use amrex_multifab_module, ONLY : amrex_mfiter, &
                                     amrex_mfiter_build, &
                                     amrex_mfiter_destroy
 
   use Grid_tile,             ONLY : Grid_tile_t
   use gr_physicalMultifabs,  ONLY : unk

   implicit none
 
   integer,           intent(IN), value :: lev
   type(c_ptr),       intent(in), value :: tags 
   real(wp),          intent(in), value :: time
   character(c_char), intent(in), value :: tagval
   character(c_char), intent(in), value :: clearval

   type(amrex_tagboxarray) :: tag
   type(amrex_mfiter)      :: mfi                                                             
   type(amrex_box)         :: bx
   type(Grid_tile_t)       :: tileDesc

   real(wp),               contiguous, pointer :: solnData(:,:,:,:)
   character(kind=c_char), contiguous, pointer :: tagData(:,:,:,:)

   integer :: refine_to
   integer :: i, j

   nullify(solnData)
   nullify(tagData)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Started on level ", lev + 1

   tag = tags

   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)                             
   do while(mfi%next())
      bx = mfi%tilebox()

      ! Level must be 1-based index and limits/limitsGC must be 1-based also
      tileDesc%level = lev + 1
      tileDesc%grid_index = mfi%grid_index()
      tileDesc%tile_index = 0
      tileDesc%limits(LOW,  :) = 1
      tileDesc%limits(HIGH, :) = 1
      tileDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
      tileDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
      tileDesc%blkLimitsGC(LOW,  :) = 1
      tileDesc%blkLimitsGC(HIGH, :) = 1
      tileDesc%blkLimitsGC(LOW,  1:NDIM) = tileDesc%limits(LOW,  1:NDIM) - NGUARD
      tileDesc%blkLimitsGC(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) + NGUARD

      call tileDesc%getDataPtr(solnData, CENTER)

      associate (lo => tileDesc%limits(LOW,  :), &
                 hi => tileDesc%limits(HIGH, :))
        tagData => tag%dataptr(mfi)
        tagData(:, :, :, :) = clearval
        do     j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)
                ! Using 0-based AMReX level indexing
                !
                ! Stored data is density.  We need to reverse transformation
                ! back to level of refinement integer that is the same on
                ! all levels
                refine_to = INT(solnData(i,j,1,1) * 4.0**(4-(lev+1)))
                if (solnData(i,j,1,1) > 0.0d0) then
                    write(*,'(A,I3,A,I3,A,F7.5,A,I3)') &
                          "     Non-zero data at (", &
                          i, ",", j, ") / density = ", solnData(i,j,1,1), &
                          " / refine to level ", refine_to
                end if

                if (refine_to > (lev+1)) then
                    write(*,'(A,I3,A,I3,A,F7.5)') "     Tag cell at (", &
                          i, ",", j, ") for refinement"
                    ! AMReX uses 0-based spatial indices/FLASH uses 1-based
                    tagData(i-1, j-1, 1, 1) = tagval
                end if
            end do
        end do
      end associate

      call tileDesc%releaseDataPtr(solnData, CENTER)
   end do
   call amrex_mfiter_destroy(mfi)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

