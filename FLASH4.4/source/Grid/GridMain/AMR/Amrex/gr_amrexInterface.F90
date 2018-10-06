!!****if* source/Grid/GridMain/AMR/Amrex/gr_amrexInterface
!!
!! NAME
!!
!!  gr_amrexInterface
!!
!! SYNOPSIS
!!
!!  use gr_amrexInterface
!!
!! DESCRIPTION 
!!
!!  Header file that defines the interfaces for all local subroutines
!!  used by AMReX-specific Grid unit code.
!!  
!!***

#include "constants.h"

module gr_amrexInterface
  interface
    subroutine gr_amrexInit()
      implicit none
    end subroutine gr_amrexInit
  end interface
  
  interface
    subroutine gr_amrexFinalize()
      implicit none
    end subroutine gr_amrexFinalize
  end interface
  
  interface
     subroutine gr_initNewLevelCallback(lev, time, pba, pdm) bind(c)
       use iso_c_binding,     ONLY : c_ptr
       use amrex_fort_module, ONLY : wp => amrex_real
       implicit none
       integer,     intent(IN), value :: lev
       real(wp),    intent(IN), value :: time
       type(c_ptr), intent(IN), value :: pba
       type(c_ptr), intent(IN), value :: pdm
     end subroutine gr_initNewLevelCallback
  end interface
  
  interface
    subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
      use iso_c_binding,     ONLY : c_ptr
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      integer,     intent(IN), value :: lev
      real(wp),    intent(IN), value :: time
      type(c_ptr), intent(IN), value :: pba
      type(c_ptr), intent(IN), value :: pdm
    end subroutine gr_makeFineLevelFromCoarseCallback
  end interface
  
  interface
    subroutine gr_remakeLevelCallback(lev, time, pba, pdm) bind(c)
      use iso_c_binding,     ONLY : c_ptr
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      integer,     intent(IN), value :: lev
      real(wp),    intent(IN), value :: time
      type(c_ptr), intent(IN), value :: pba
      type(c_ptr), intent(IN), value :: pdm
    end subroutine gr_remakeLevelCallback
  end interface 
  
  interface
    subroutine gr_clearLevelCallback(lev) bind(c)
      implicit none
      integer, intent(IN), value :: lev
    end subroutine gr_clearLevelCallback
  end interface
  
  interface
    subroutine gr_markRefineDerefineCallback(lev, tags, time, &
                                           tagval, clearval) bind(c)
      use iso_c_binding,     ONLY : c_ptr, c_char
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      integer,           intent(IN), value :: lev
      type(c_ptr),       intent(IN), value :: tags 
      real(wp),          intent(IN), value :: time
      character(c_char), intent(IN), value :: tagval
      character(c_char), intent(IN), value :: clearval
    end subroutine gr_markRefineDerefineCallback
  end interface
 
  interface
    subroutine gr_fillPhysicalBC(pmf, scomp, ncomp, time, pgeom) bind(c)
      use iso_c_binding,     ONLY : c_ptr, c_int
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      type(c_ptr),    value :: pmf
      type(c_ptr),    value :: pgeom
      integer(c_int), value :: scomp
      integer(c_int), value :: ncomp
      real(wp),       value :: time
    end subroutine gr_fillPhysicalBC
  end interface

  interface
    subroutine gr_fillPhysicalFaceBC(pmf, scomp, ncomp, time, pgeom) bind(c)
      use iso_c_binding,     ONLY : c_ptr, c_int
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      type(c_ptr),    value :: pmf
      type(c_ptr),    value :: pgeom
      integer(c_int), value :: scomp
      integer(c_int), value :: ncomp
      real(wp),       value :: time
    end subroutine gr_fillPhysicalFaceBC
  end interface

  interface
    subroutine gr_getFinestLevel(level)
      implicit none
      integer, intent(IN) :: level
    end subroutine gr_getFinestLevel
  end interface

  interface
    subroutine gr_writeData(stepno, t_new, argBaseName)
      implicit none
      integer, intent(IN) :: stepno
      real,    intent(IN) :: t_new
      character(len=*), intent(IN), optional :: argBaseName
    end subroutine gr_writeData
  end interface
 
  interface
    subroutine gr_preinterpolationWork(lo, hi, &
                                       d, dlo, dhi, nd, &
                                       scomp, ncomp) bind(c)
      use iso_c_binding,     ONLY : c_int
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      integer(c_int), intent(in)          :: lo(MDIM), hi(MDIM)
      integer(c_int), intent(in)          :: dlo(MDIM), dhi(MDIM)
      integer(c_int), intent(in),   value :: nd
      integer(c_int), intent(in),   value :: scomp
      integer(c_int), intent(in),   value :: ncomp
      real(wp),       intent(inout)       :: d(dlo(IAXIS):dhi(IAXIS), &
                                               dlo(JAXIS):dhi(JAXIS), &
                                               dlo(KAXIS):dhi(KAXIS), &
                                               nd)
    end subroutine gr_preinterpolationWork
  end interface

  interface
    subroutine gr_postinterpolationWork(lo, hi, &
                                        d, dlo, dhi, nd, &
                                        scomp, ncomp) bind(c)
      use iso_c_binding,     ONLY : c_int
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      integer(c_int), intent(in)          :: lo(MDIM), hi(MDIM)
      integer(c_int), intent(in)          :: dlo(MDIM), dhi(MDIM)
      integer(c_int), intent(in),   value :: nd
      integer(c_int), intent(in),   value :: scomp
      integer(c_int), intent(in),   value :: ncomp
      real(wp),       intent(inout)       :: d(dlo(IAXIS):dhi(IAXIS), &
                                               dlo(JAXIS):dhi(JAXIS), &
                                               dlo(KAXIS):dhi(KAXIS), &
                                               nd)
    end subroutine gr_postinterpolationWork
  end interface

  interface
    subroutine gr_primitiveToConserve(lo, hi, &
                                      d, dlo, dhi, nd, &
                                      scomp, ncomp)
      implicit none
      integer, intent(in)    :: lo(MDIM), hi(MDIM)
      integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
      integer, intent(in)    :: nd
      integer, intent(in)    :: scomp
      integer, intent(in)    :: ncomp
      real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                                  dlo(JAXIS):dhi(JAXIS), &
                                  dlo(KAXIS):dhi(KAXIS), &
                                  nd)
    end subroutine gr_primitiveToConserve
  end interface

  interface
    subroutine gr_conserveToPrimitive(lo, hi, &
                                      d, dlo, dhi, nd, &
                                      scomp, ncomp)
      implicit none
      integer, intent(in)    :: lo(MDIM), hi(MDIM)
      integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
      integer, intent(in)    :: nd
      integer, intent(in)    :: scomp
      integer, intent(in)    :: ncomp
      real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                                  dlo(JAXIS):dhi(JAXIS), &
                                  dlo(KAXIS):dhi(KAXIS), &
                                  nd)
    end subroutine gr_conserveToPrimitive
  end interface
 
  interface
    subroutine gr_cleanDensityData(smallRho, &
                                   lo, hi, &
                                   d, dlo, dhi, nd)
      implicit none
      real,    intent(in)    :: smallRho
      integer, intent(in)    :: lo(MDIM), hi(MDIM)
      integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
      integer, intent(in)    :: nd
      real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                                  dlo(JAXIS):dhi(JAXIS), &
                                  dlo(KAXIS):dhi(KAXIS), &
                                  nd)
    end subroutine gr_cleanDensityData
  end interface

  interface
    subroutine gr_cleanEnergyData(smallE, &
                                   lo, hi, &
                                   d, dlo, dhi, nd)
      implicit none
      real,    intent(in)    :: smallE
      integer, intent(in)    :: lo(MDIM), hi(MDIM)
      integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
      integer, intent(in)    :: nd
      real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                                  dlo(JAXIS):dhi(JAXIS), &
                                  dlo(KAXIS):dhi(KAXIS), &
                                  nd)
    end subroutine gr_cleanEnergyData
  end interface

  interface
    subroutine gr_restrictAllLevels(gridDataStruct, convertPtoC, convertCtoP)
      implicit none
      integer, intent(IN) :: gridDataStruct
      logical, intent(IN) :: convertPtoC
      logical, intent(IN) :: convertCtoP
    end subroutine gr_restrictAllLevels
  end interface
 
  interface
    subroutine gr_splitFabAtBoundary(face, axis, patch, delta, &
                                     interior, guardcells, hasGC)
      implicit none
      integer, intent(IN)  :: face
      integer, intent(IN)  :: axis
      integer, intent(IN)  :: patch(LOW:HIGH, 1:MDIM)
      real,    intent(IN)  :: delta(1:MDIM)
      integer, intent(OUT) :: interior(LOW:HIGH, 1:MDIM)
      integer, intent(OUT) :: guardcells(LOW:HIGH, 1:MDIM)
      logical, intent(OUT) :: hasGC
    end subroutine gr_splitFabAtBoundary
  end interface

  interface
    subroutine gr_copyFabInteriorToRegion(fab, gds, face, axis, interior, &
                                          scomp, ncomp, region)
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      real(wp), pointer, contiguous, intent(IN)    :: fab(:, :, :, :)
      integer,                       intent(IN)    :: gds
      integer,                       intent(IN)    :: face
      integer,                       intent(IN)    :: axis
      integer,                       intent(IN)    :: interior(LOW:HIGH, 1:MDIM)
      integer,                       intent(IN)    :: scomp
      integer,                       intent(IN)    :: ncomp
      real(wp), pointer, contiguous, intent(INOUT) :: region(:, :, :, :)
    end subroutine gr_copyFabInteriorToRegion
  end interface
 
  interface
    subroutine gr_copyGuardcellRegionToFab(region, gds, face, axis, guardcells, &
                                           scomp, ncomp, fab)
      use amrex_fort_module, ONLY : wp => amrex_real
      implicit none
      real(wp), intent(IN),    pointer, contiguous :: region(:, :, :, :)
      integer,  intent(IN)                         :: gds
      integer,  intent(IN)                         :: face
      integer,  intent(IN)                         :: axis
      integer,  intent(IN)                         :: guardcells(LOW:HIGH, 1:MDIM)
      integer,                       intent(IN)    :: scomp
      integer,                       intent(IN)    :: ncomp
      real(wp), intent(INOUT), pointer, contiguous :: fab(:, :, :, :)
    end subroutine gr_copyGuardcellRegionToFab
  end interface

end module gr_amrexInterface

