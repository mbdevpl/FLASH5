!!****if* source/physics/IncompNS/IncompNSMain/vardens/ins_computeDtLocal
!!
!! NAME
!!
!!  ins_computeDtLocal
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***

subroutine ins_computeDtLocal(blockID,   & 
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc )

  use IncompNS_data, ONLY : ins_cflflg, ins_cfl, ins_sigma, ins_invRe, ins_dtspec

  use Grid_interface, ONLY : Grid_getBlkCenterCoords
  
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash.h"
#include "constants.h"
  integer, intent(IN) :: blockID
  integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  integer, intent(IN) :: isize,jsize,ksize
  real, intent(IN) :: dx, dy, dz
  real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  ! Local variables:
  real, parameter :: eps = 1.e-12
  real :: dtc,dtv,dtl,velcoeff


  if (ins_cflflg .eq. 0) then
     dtlocal    = ins_dtspec
     lminloc(:) = 0
     return
  endif

# if NDIM == MDIM

  velcoeff =  MAX( MAXVAL(ABS(facexData(VELC_FACE_VAR,:,:,:))/dx), &
                   MAXVAL(ABS(faceyData(VELC_FACE_VAR,:,:,:))/dy), &
                   MAXVAL(ABS(facezData(VELC_FACE_VAR,:,:,:))/dz) )

  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps
  endif
  
  dtv = ins_sigma / (ins_invRe*MAX( 1.0/(dx*dx), 1.0/(dy*dy),1.0/(dz*dz) ))
         
# elif NDIM == 2

  velcoeff =  MAX( MAXVAL(ABS(facexData(VELC_FACE_VAR,:,:,:))/dx), &
                   MAXVAL(ABS(faceyData(VELC_FACE_VAR,:,:,:))/dy))

  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps  
  endif
  
  dtv = ins_sigma / (ins_invRe*MAX( 1.0/(dx*dx), 1.0/(dy*dy)))

# endif

  dtl = 0.5*MIN(dtc,dtv)

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
     ! Cell located at center of Block - Used to define block location. 
     lminloc(1) = NGUARD + NXB/2
     lminloc(2) = NGUARD + NYB/2
     lminloc(3) = NGUARD + NZB/2
     lminloc(4) = blockID
     lminloc(5) = gr_meshMe

  endif

end subroutine ins_computeDtLocal
