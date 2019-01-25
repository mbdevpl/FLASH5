!!****if* source/physics/Hydro/HydroMain/unsplit/hy_addGravity
!!
!!
!! NAME
!!
!!  hy_addGravity
!!
!! SYNOPSIS
!!
!!  call hy_addGravity( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimits(LOW:HIGH,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real(IN)    :: dt,
!!                            real(IN)    :: gravX,
!!                            real(IN)    :: gravY,
!!                            real(IN)    :: gravZ )
!!
!!
!! DESCRIPTION
!!
!!  Adds the second part of the gravitational force to the momenta and energy.
!!  The first half is added in unsplitUpdate.  This centers the gravitational
!!  force appropriately on the current timestep rather than extrapolating.
!!
!! ARGUMENTS
!!
!!  blockID   - local block ID
!!  blkLimits - block limits 
!!  dataSize  - array size of gravX,Y,Z in non-fixed block size mode in UG
!!  dt        - timestep
!!  gravX     - gravity source term in x-direction
!!  gravY     - gravity source term in y-direction
!!  gravZ     - gravity source term in z-direction
!!
!! NOTES
!!
!!  The ENER_VAR component of UNK is now always in specific
!!  (i.e., energy per mass) form, on entry as well as on return from this routine.
!!  This is changed from the behavior in FLASH4.2.2 and earlier.
!!
!!***

!!REORDER(4):U

Subroutine hy_addGravity&
     (tileDesc,blkLimits,loGC,hiGC,dt,gravX,gravY,gravZ)

  use Hydro_data,      ONLY : hy_useGravity,        &
                              hy_useGravHalfUpdate

  use Grid_interface,  ONLY : Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr
  use flash_tile, ONLY : flash_tile_t

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  type(flash_tile_t), intent(IN)   :: tileDesc

  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  integer, intent(IN), dimension(MDIM):: loGC, hiGC
  real,    intent(IN) :: dt

  real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), intent(IN) :: &
       gravX,gravY,gravZ
  !! -----------------------------------------------------
  real :: hdt
  integer :: i,j,k
  real, pointer, dimension(:,:,:,:) :: U
  real, dimension(3) :: velNew

  hdt = 0.5 * dt

  !! Get block pointer for storages of Riemann states
  call tileDesc%getDataPtr(U, CENTER)

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) .LE. 0.0) then
#endif
              velNew(1:3) = U(VELX_VAR:VELZ_VAR,i,j,k)&
                + hdt*(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/)

              U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) &
                + hdt*dot_product(velNew(1:3),(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/))

              U(VELX_VAR:VELZ_VAR,i,j,k) = velNew(1:3)
#ifdef BDRY_VAR
           endif
#endif
        enddo
     enddo
  enddo

  !! Release block pointer for storages of Riemann states
  call tileDesc%releaseDataPtr(U, CENTER)

end Subroutine hy_addGravity
