subroutine gr_makeNewLevelFromScratch(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,      ONLY : wp => amrex_real
    use amrex_amr_module,       ONLY : amrex_geom, &
                                       amrex_problo
    use amrex_box_module,       ONLY : amrex_box
    use amrex_boxarray_module,  ONLY : amrex_boxarray, &
                                       box_print => amrex_print
    use amrex_distromap_module, ONLY : amrex_distromap, &
                                       distro_print => amrex_print
    use amrex_parallel_module,  ONLY : amrex_parallel_myproc
    use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                       amrex_mfiter_build, &
                                       amrex_mfiter_destroy, &
                                       amrex_multifab_build
    use physicaldata,           ONLY : unk

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

    rank = amrex_parallel_myproc()
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - Start Level ", lev
 
    ba = pba
    dm = pdm
 
    call box_print(ba)
    call distro_print(dm)

    call gr_clearLevel(lev)

    ! Create FABS for storing physical data at coarsest level
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)

    ! Write initial data across domain at coarsest level
    call amrex_mfiter_build(mfi, unk(lev))

    do while (mfi%next())
        bx = mfi%tilebox()
        initData => unk(lev)%dataptr(mfi)

        associate(lo => bx%lo, &
                  hi => bx%hi, &
                  dx => amrex_geom(lev)%dx)
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
 
    write(*,*) "[Rank ", rank, "] gr_makeNewLevelFromScratch - End Level ", lev
end subroutine gr_makeNewLevelFromScratch 

