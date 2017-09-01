subroutine gr_makeNewLevelFromScratch(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom, &
                                       amrex_problo
    use amrex_box_module,       ONLY : amrex_box
    use amrex_boxarray_module,  ONLY : amrex_boxarray
    use amrex_distromap_module, ONLY : amrex_distromap
    use amrex_parallel_module,  ONLY : amrex_parallel_myproc
    use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy, &
                                       amrex_multifab_build
    
    use gr_physicalMultifabs,   ONLY : unk, &
                                       facevarx, facevary, facevarz
    use amrex_interfaces,       ONLY : gr_clearLevel

    implicit none

#include "constants.h"
#include "Flash.h"

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter)    :: mfi
    type(amrex_box)       :: bx

    real(wp), contiguous, pointer :: initData(:,:,:,:)
    real(wp)                      :: x = 0.0_wp 
    real(wp)                      :: y = 0.0_wp
    real(wp)                      :: z = 0.0_wp

    integer :: rank = 0

    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    
    integer :: lev_flash = -1

    ! C++ AMReX uses zero-based level index set, but FLASH uses 1-based set
    lev_flash = lev + 1 

    rank = amrex_parallel_myproc()
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - Start Level ", &
                lev_flash
 
    ba = pba
    dm = pdm

    call gr_clearLevel(lev)

    ! Create FABS for storing physical data at coarsest level
    call amrex_multifab_build(unk(lev_flash), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarx(lev_flash), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevary(lev_flash), ba, dm, NUNK_VARS, NGUARD)
    call amrex_multifab_build(facevarz(lev_flash), ba, dm, NUNK_VARS, NGUARD)

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev_flash))

    do while (mfi%next())
        bx = mfi%tilebox()
        initData => unk(lev_flash)%dataptr(mfi)

        associate(lo => bx%lo, &
                  hi => bx%hi, &
                  dx => amrex_geom(lev_flash)%dx)
            do         k = lo(3), hi(3) 
                do     j = lo(2), hi(2)
                    z = amrex_problo(3) + (dble(k)+0.5d0) * dx(3)
                    y = amrex_problo(2) + (dble(j)+0.5d0) * dx(2)
                    do i = lo(1), hi(1)
                        x = amrex_problo(1) + (dble(i)+0.5d0) * dx(1)
                        initData(i, j, k, UNK_VARS_BEGIN) = x + y + z 
                    end do
                end do
            end do
        end associate
    end do

    call amrex_mfiter_destroy(mfi)
 
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - End Level ", &
                lev_flash
end subroutine gr_makeNewLevelFromScratch 

