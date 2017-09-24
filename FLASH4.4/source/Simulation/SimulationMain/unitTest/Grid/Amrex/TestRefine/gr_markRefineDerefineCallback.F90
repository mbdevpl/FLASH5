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

        tagData(:, :, :, :) = clearval
        do     j = lo(2), hi(2)
            do i = lo(1), hi(1)
                ! Using 0-based AMReX level indexing
                !
                ! Stored data is density.  We need to reverse transformation
                ! back to level of refinement integer that is the same on
                ! all levels
                refine_to = INT(solnData(i,j,1,1) * 4.0d0**(4-(lev+1)))
                if (solnData(i,j,1,1) > 0.0d0) then
                    write(*,'(A,I2,A,I2,A,F7.5,A,I2)') &
                          "     Non-zero data at (", &
                          i, ",", j, ") / density = ", solnData(i,j,1,1), &
                          " / refine to level ", refine_to
                end if

                if (refine_to > (lev+1)) then
                    write(*,'(A,I2,A,I2,A,F7.5)') "     Tag cell at (", &
                          i, ",", j, ") for refinement"
                    tagData(i, j, 1, 1) = tagval
                end if
            end do
        end do
      end associate
   end do
   call amrex_mfiter_destroy(mfi)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

