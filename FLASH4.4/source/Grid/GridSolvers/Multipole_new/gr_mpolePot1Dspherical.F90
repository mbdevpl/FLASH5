!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePot1Dspherical
!!
!! NAME
!!
!!  gr_mpolePot1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot1Dspherical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a one-dimensional spherical geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar : index to variable containing the potential
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot1Dspherical (ipotvar)

  use Grid_interface,    ONLY : Grid_getBlkPtr,        &
                                Grid_releaseBlkPtr,    &
                                Grid_getBlkBoundBox,   &
                                Grid_getDeltas,        &
                                Grid_getLeafIterator,   &
                                Grid_releaseLeafIterator

  use gr_mpoleData,      ONLY : gr_mpoleGravityConstant,        &
                                gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxRadialZones,         & 
                                gr_mpoleMinRadialZone,          & 
                                gr_mpoleZoneRmax,               &
                                gr_mpoleZoneQmax,               &
                                gr_mpoleZoneType,               &
                                gr_mpoleZoneScalarInv,          &
                                gr_mpoleZoneLogNormInv,         &
                                gr_mpoleZoneExponentInv,        &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneQshift,        &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI

  use block_metadata,    ONLY : block_metadata_t
  use leaf_iterator,     ONLY : leaf_iterator_t

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: ipotvar

  logical :: innerZonePotential

  
  integer :: DrUnit
  integer :: i, n
  integer :: imax, imin
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: zone

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  

  real    :: bndBoxILow
  real    :: DeltaI, DeltaIHalf, DeltaIFourth
  real    :: potential
  real    :: Qfloat, QfracI, QfracR
  real    :: Rcenter, Rsph
  real    :: RdotI, IdotR
  real    :: rlocal, rinDrs, rinvI
  real    :: sclInv, lgnInv, expInv

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, pointer :: solnData (:,:,:,:)
!  
  integer :: lev
  type(block_metadata_t) :: block
  type(leaf_iterator_t) :: itor
!
!     ...Sum quantities over all locally held leaf blocks.
!
!

  ! Replaced `!$omp do schedule(static)` with `!$omp single` below as temporary fix until we determine
  ! the proper way to parallelize leaf iterator loops with OpenMP - JAH
  !$omp single
  call Grid_getLeafIterator(itor)
  do while(itor%is_valid())
     call itor%blkMetaData(block)
     lev=block%level
     blkLimits=block%limits
     
     call Grid_getBlkBoundBox     (block, bndBox)
     call Grid_getDeltas          (lev, delta)
     call Grid_getBlkPtr          (block, solnData)

     imin = blkLimits (LOW, IAXIS)
     imax = blkLimits (HIGH,IAXIS)

     DeltaI       = delta (IAXIS)
     DeltaIHalf   = DeltaI * HALF
     DeltaIFourth = DeltaI * FOURTH

     bndBoxILow = bndBox (LOW,IAXIS)

     solnData (ipotvar , imin:imax , 1,1) = ZERO
!
!
!          ...The 1D spherical case:
!
!
!                                    |                   |            O --> multipole expansion origin
!                O----n----i----n----|----n----i----n----|---         i --> cell center
!                                    |                   |            n --> Rsph: 1/4 cell size off center
!
!
!             The potentials will not be evaluated at the cell centers but rather
!             at 2 points off 1/4 cell size off the center's position. It is not
!             possible to use the cell faces for potential averaging, due to the 
!             r = 0 condition at the multipole expansion origin. Only L = 0
!             solid harmonics are calculated and combined with the moments.
!
!
     Rcenter = bndBoxILow + DeltaIHalf

     do i = imin, imax                                  ! loop over all cells in block

        do n = -1,1,2                                   ! loop over the two points off the cell center

           Rsph = Rcenter + real (n) * DeltaIFourth     ! radius off by 1/4 cell size from cell center
!
!
!        ...Find the radial bin.
!
!
           innerZonePotential = Rsph <= gr_mpoleInnerZoneMaxR


           if (innerZonePotential) then

               rinDrs = Rsph * gr_mpoleDrInnerZoneInv
               DrUnit = int (ceiling (rinDrs))
               Qlower = gr_mpoleInnerZoneQlower (DrUnit)
               Qupper = gr_mpoleInnerZoneQupper (DrUnit)
               QfracR = ZERO
               QfracI = ONE

               do Q = Qlower,Qupper
                  if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
               end do

           else

               do zone = gr_mpoleMinRadialZone , gr_mpoleMaxRadialZones
                  if (Rsph - gr_mpoleZoneRmax (zone) <= ZERO) exit
               end do

               rlocal = Rsph - gr_mpoleZoneRmax (zone - 1)
               type   = gr_mpoleZoneType        (zone)
               sclInv = gr_mpoleZoneScalarInv   (zone)
               expInv = gr_mpoleZoneExponentInv (zone)

               if (type == ZONE_EXPONENTIAL) then
                   Qfloat = (rlocal * sclInv * gr_mpoleDrInv) ** expInv
               else if (type == ZONE_LOGARITHMIC) then
                   lgnInv = gr_mpoleZoneLogNormInv (zone)
                   Qfloat = expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE)
               end if

               Qlocal = ceiling (Qfloat)
               QfracI = real (Qlocal) - Qfloat
               QfracR = ONE - QfracI
               Q      = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

           end if
!
!
!        ...Calculate and add the current potential to the current cell.
!
!
           rinvI = ONE / Rsph

           RdotI =   QfracI * gr_mpoleQDampingI (Q)   * gr_mpoleMomentI (1,Q  )   &
                   + QfracR * gr_mpoleQDampingI (Q+1) * gr_mpoleMomentI (1,Q+1)

           IdotR = rinvI * (  QfracI * gr_mpoleMomentR (1,Q-1) &
                            + QfracR * gr_mpoleMomentR (1,Q)   )

           potential = - gr_mpoleGravityConstant * (RdotI + IdotR)

           solnData (ipotvar,i,1,1) = solnData (ipotvar,i,1,1) + potential

        end do
        Rcenter = Rcenter + DeltaI
     end do
!
!
!    ...Form the potential average in each cell.
!
!
     solnData (ipotvar,imin:imax,1,1) = HALF * solnData (ipotvar,imin:imax,1,1)
!
!
!    ...Get ready for retrieving next LEAF block for the current processor.
!
!
     call Grid_releaseBlkPtr (block, solnData)
     call itor%next()
  end do
  call Grid_releaseLeafIterator(itor)
  !$omp end single
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpolePot1Dspherical

