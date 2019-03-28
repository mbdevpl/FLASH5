!!****f* source/Grid/GridMain/AMR/Amrex/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!  call Grid_conserveFluxes(integer(IN) :: axis,
!!                           integer(IN) :: coarse_level)
!!  
!! DESCRIPTION 
!!  When FLASH is run with AMR, it is possible that some leaf blocks
!!  will have a neighboring leaf block that is refined at the next coarsest
!!  level.  To maintain conservation, the flux entering into the
!!  fine block at this shared boundary must equal the flux leaving the
!!  coarse block at the boundary.
!!
!!  To enforce this requirement, this routine overwrites in the coarse
!!  block the flux at the shared boundary with the averaged flux from the
!!  fine block, which is sensible as the fine flux data should be at least as
!!  accurate as the flux in the coarse level.
!! 
!!  It is assumed that before calling this routine, the code has already
!!  loaded the corrected fluxes for the fine level into the flux 
!!  registers using Grid_putFluxData and that the uncorrected flux for the
!!  coarse level has been stored in Grid with the Grid_getFluxPtr interface.
!!  This routine will only overwrite the flux data at fine/coarse boundaries.
!!
!! ARGUMENTS 
!!  axis - the only acceptable value for AMReX is ALLDIR.
!!  coarse_level - the 1-based level index of the coarse blocks
!!
!! SEE ALSO
!!  Grid_getFluxPtr/Grid_releaseFluxPtr
!!  Grid_putFluxData
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"
#include "constants.h"

subroutine Grid_conserveFluxes(axis, coarse_level, isDensity)
    use amrex_fort_module,    ONLY : wp => amrex_real
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    use Driver_interface,     ONLY : Driver_abortFlash
    use Grid_interface,       ONLY : Grid_getGeometry, &
                                     Grid_getTileIterator, &
                                     Grid_releaseTileIterator, &
                                     Grid_getCellFaceAreas
    use Grid_iterator,        ONLY : Grid_iterator_t
    use Grid_tile,            ONLY : Grid_tile_t
    use gr_physicalMultifabs, ONLY : fluxes, &
                                     flux_registers

    implicit none

    integer, intent(IN)                   :: axis
    integer, intent(IN)                   :: coarse_level
    integer, intent(IN), optional, target :: isDensity(:)

    integer :: fine
    integer :: coarse
    integer :: geometry

    real, pointer     :: fluxData(:,:,:,:)
    real, allocatable :: faceAreas(:,:,:)

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t)     :: tileDesc

    integer :: lo(4)
    integer :: hi(4)

    integer :: i, j, k, var

    nullify(fluxData)

    if (axis /= ALLDIR) then
        call Driver_abortFlash("[Grid_conserveFluxes] AMReX requires axis==ALLDIR")
    end if
    ! DEV: Accept this for now as Hydro is calling it, but ignore the contents
!    if (present(isDensity)) then
!        call Driver_abortFlash("[Grid_conserveFluxes] isDensity not implemented")
!    end if
    
    ! FLASH uses 1-based level index / AMReX uses 0-based index
    coarse = coarse_level - 1
    fine   = coarse_level

    ! No need to conserve on the finest level in existence or any
    ! level index corresponding to a finer mesh
    if (coarse >= amrex_get_finest_level())     RETURN

    call Grid_getTileIterator(itor, ALL_BLKS, level=coarse_level, tiling=.FALSE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       call tileDesc%getDataPtr(fluxData, FLUXX)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(IAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                   ! Basic sanity check of faceAreas
                   if (geometry == CYLINDRICAL) then
                      if ((i == 1) .AND. (faceAreas(i,j,k) /= 0.0)) then
                         write(*,*) "face area != 0 for r==0", i, j, k
                         STOP
                      end if
                      if ((faceAreas(i,j,k) == 0.0) .AND. (i /= 1)) then
                         write(*,*) "Zero face area for r > 0", i, j, k
                         STOP
                      end if
                   end if
#endif

                   ! There is potentially non-zero flux density data at r=0 that
                   ! should remain unchanged during the whole flux correction
                   ! process.  Therefore, don't transform the r=0 flux density data.  
                   if (faceAreas(i, j, k) /= 0.0) then
                      fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                               * faceAreas(i, j, k)
                   end if
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXX)

#if   NDIM >= 2
       call tileDesc%getDataPtr(fluxData, FLUXY)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(JAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                   if (faceAreas(i,j,k) == 0.0) then
                      write(*,*) "Zero face area along J at", i, j, k
                      STOP
                   end if
#endif
                   fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                            * faceAreas(i, j, k)
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXY)
#endif

#if   NDIM == 3
       call tileDesc%getDataPtr(fluxData, FLUXZ)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(KAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                    if (faceAreas(i,j,k) == 0.0) then
                       write(*,*) "Zero face area along K at", i, j, k
                       STOP
                    end if
#endif
                    fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                             * faceAreas(i, j, k)
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXZ)
#endif

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    ! This should not overwrite any of the r=0 flux density data as those
    ! faces are at the domain boundary rather than a fine/coarse boundary
    call flux_registers(fine)%overwrite(fluxes(coarse, :), 1.0_wp)

    call Grid_getTileIterator(itor, ALL_BLKS, level=coarse_level, tiling=.FALSE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       call tileDesc%getDataPtr(fluxData, FLUXX)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(IAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                   ! Basic sanity check of faceAreas
                   if (geometry == CYLINDRICAL) then
                      if ((i == 1) .AND. (faceAreas(i,j,k) /= 0.0)) then
                         write(*,*) "face area != 0 for r==0", i, j, k
                         STOP
                      end if
                      if ((faceAreas(i,j,k) == 0.0) .AND. (i /= 1)) then
                         write(*,*) "Zero face area for r > 0", i, j, k
                         STOP
                      end if
                   end if
#endif

                   ! Again, leave the flux density data at r=0 alone so
                   ! that the original data is untouched by the whole
                   ! flux conservation process
                   if (faceAreas(i, j, k) /= 0.0) then
                      fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                               / faceAreas(i, j, k)
                   end if
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXX)

#if   NDIM >= 2
       call tileDesc%getDataPtr(fluxData, FLUXY)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(JAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                   if (faceAreas(i,j,k) == 0.0) then
                      write(*,*) "Zero face area along J at", i, j, k
                      STOP
                   end if
#endif
                   fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                            / faceAreas(i, j, k)
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXY)
#endif

#if   NDIM == 3
       call tileDesc%getDataPtr(fluxData, FLUXZ)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       allocate(faceAreas(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(KAXIS, tileDesc%level, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

       do        var = 1, NFLUXES
          do       k = lo(KAXIS), hi(KAXIS)
             do    j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
#ifdef DEBUG_GRID
                   if (faceAreas(i,j,k) == 0.0) then
                      write(*,*) "Zero face area along K at", i, j, k
                      STOP
                   end if
#endif
                   fluxData(i, j, k, var) =    fluxData(i, j, k, var) &
                                            / faceAreas(i, j, k)
                end do
             end do
          end do
       end do

       deallocate(faceAreas)
       call tileDesc%releaseDataPtr(fluxData, FLUXZ)
#endif

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

!    call Grid_getGeometry(geometry)
!
!    select case (geometry)
!    case (CARTESIAN)
!        ! The AMReX flux registers are dealing with fluxes and 
!        ! *not* flux densities.  In Grid_putFluxData, flux densities at the fine
!        ! level were scaled to fluxes with the assumption that the cell lengths
!        ! at the coarse level are one.  Therefore, reconversion to flux densities
!        ! is automatic here.
!        call flux_registers(fine)%overwrite(fluxes(coarse, :), 1.0_wp)
!    case default
!        ! DEV: TODO This routine should take an isFluxDensity array as an argument
!        ! so that the routine can determine which need to be scaled and which do not
!        call Driver_abortFlash("[Grid_conserveFluxes] Only works with Cartesian")
!    end select
end subroutine Grid_conserveFluxes

