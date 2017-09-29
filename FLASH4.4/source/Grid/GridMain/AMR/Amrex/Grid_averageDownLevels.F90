#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"

subroutine Grid_averageDownLevels()
    use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                          amrex_geom, &
                                          amrex_ref_ratio
    use amrex_multifabutil_module, ONLY : amrex_average_down

    use gr_physicalMultifabs,      ONLY : unk

    implicit none

    integer :: lev
    integer :: finest_level

    finest_level = amrex_get_finest_level()

    ! DEV: TODO Implement for face variables as well
    do lev = finest_level, 1, -1
        call amrex_average_down(unk(lev  ), &
                                unk(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                UNK_VARS_BEGIN, NUNK_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 

#ifdef DEBUG_GRID
    if (finest_level == 0) then
        write(*,'(A,A)') "[gr_averageDownLevels]", &
                            "               No need to average"
    else
        write(*,'(A,A,I2,A)') "[gr_averageDownLevels]", &
                              "               From ", &
                             finest_level, " down to 1"
    end if
#endif

end subroutine Grid_averageDownLevels

