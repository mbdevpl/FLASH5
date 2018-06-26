!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_assembleFluidForce_toPoint
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

subroutine sm_surf_assembleFluidForce_toPoint(ibd, point, force_pres, force_visc, moment_pres, moment_visc)
  use gr_sbData, ONLY : gr_sbBodyInfo
  use sm_Misc_interface, only: sm_crossProd
  use Driver_interface, only: Driver_abortFlash
  implicit none

  ! IO variables 
  integer, intent(in) :: ibd
  real, dimension(NDIM), intent(in) :: point  ! location to compute moments about
  real, dimension(NDIM), intent(out):: force_pres, force_visc, moment_pres, moment_visc

  ! Local Variables
  integer :: p
  real, dimension(NPART_PROPS) :: particle
  real, dimension(NDIM) :: nhat, f_pres, f_visc, pos, rpcm
  real                  :: area
  
  ! Loop over particles
  force_pres  = 0.
  force_visc  = 0.
  moment_pres = 0.
  moment_visc = 0.
  do p = 1,gr_sbBodyInfo(ibd)%totalPart
     
     ! get patch info
     particle = gr_sbBodyInfo(ibd)% particles(1:NPART_PROPS,p)
     
     ! get the unit norm of the patch
#if NDIM == 2
     nhat = (/ particle(NMLX_PART_PROP), particle(NMLY_PART_PROP) /)
#elif NDIM == 3
     nhat = (/ particle(NMLX_PART_PROP), particle(NMLY_PART_PROP), particle(NMLZ_PART_PROP) /)
#endif

     ! get the force on the patch from the particle
     ! $ f = \vec{f}_\nu - pres*\hat{n} $
     ! Pressure traction
     f_pres = - particle(PRES_PART_PROP)*nhat
     ! Visc Traction
#if NDIM == 2
     f_visc = (/ particle(FXVI_PART_PROP), particle(FYVI_PART_PROP) /) 
#elif NDIM == 3
     f_visc = (/ particle(FXVI_PART_PROP), particle(FYVI_PART_PROP), particle(FZVI_PART_PROP) /)
#endif
     
     ! area
     area = particle(AREA_PART_PROP)

     ! position of center of patch
#if NDIM == 2
     pos = (/ particle(POSX_PART_PROP), particle(POSY_PART_PROP) /)
#elif NDIM == 3
     pos = (/ particle(POSX_PART_PROP), particle(POSY_PART_PROP), particle(POSZ_PART_PROP) /)
#endif

     ! compute total force
     force_pres = force_pres + f_pres*area
     force_visc = force_visc + f_visc*area

     rpcm(1:NDIM) = pos(1:NDIM)-point(1:NDIM)
     ! compute total moment
#if NDIM == 2
     ! M = M + rx*Fy - ry*Fx
     !moment_pres(CONSTANT_ONE) = moment_pres(CONSTANT_ONE)     + &
     !                           (rpcm(IAXIS)*force_pres(JAXIS) - &
     !                            rpcm(JAXIS)*force_pres(IAXIS))*area   
     !moment_visc(CONSTANT_ONE) = moment_visc(CONSTANT_ONE)     + &
     !                           (rpcm(IAXIS)*force_visc(JAXIS) - &
     !                            rpcm(JAXIS)*force_visc(IAXIS))*area
 
     moment_pres(CONSTANT_ONE) = moment_pres(CONSTANT_ONE)     + &
                                 sm_crossProd(rpcm, f_pres*area )
     moment_visc(CONSTANT_ONE) = moment_visc(CONSTANT_ONE)     + &
                                 sm_crossProd(rpcm, f_visc*area )

#elif NDIM == 3
     moment_pres = moment_pres + sm_crossProd(rpcm, f_pres*area ) 
     moment_visc = moment_visc + sm_crossProd(rpcm, f_visc*area ) 
#endif   
  
  end do

  return

end subroutine sm_surf_assembleFluidForce_toPoint
