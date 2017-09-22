subroutine gr_averageDownLevels()
    use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                          amrex_geom, &
                                          amrex_ref_ratio
    use amrex_multifabutil_module, ONLY : amrex_average_down

    use gr_physicalMultifabs,      ONLY : unk

    implicit none

#include "Flash.h"

    integer :: lev
    integer :: finest_level

    finest_level = amrex_get_finest_level()

    ! DEV: TODO Implement for face variables as well
    do lev = finest_level-1, 0, -1
        call amrex_average_down(unk(lev+1), &
                                unk(lev  ), &
                                amrex_geom(lev+1), &
                                amrex_geom(lev  ), &
                                UNK_VARS_BEGIN, NUNK_VARS, &
                                amrex_ref_ratio(lev))
    end do 

    write(*,'(A,I2,A)') "[gr_averageDownLevels]          From ", &
                        finest_level, " down to 1"
end subroutine gr_averageDownLevels

