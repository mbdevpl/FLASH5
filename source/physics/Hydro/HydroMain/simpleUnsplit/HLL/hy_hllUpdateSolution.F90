! Note: the following arrays need to be spelled exactly like this in the code below,
!       preserving case.
!!REORDER(4): Uin, Uout, fl[XYZ]

#ifdef DEBUG_ALL
#define DEBUG_UHD
#endif

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
#include "UHD.h"

Subroutine hy_hllUpdateSolution( tileLimits, Uin, plo, Uout, flX, flY, flZ, loFl, del, dt )
  use Hydro_data,        ONLY : hy_fluxCorrect,      &
                                hy_useGravity,       &
                                hy_unsplitEosMode,   &
                                hy_updateHydroFluxes
  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampVarMask

  implicit none

  !! ---- Argument List ----------------------------------
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: UOUT(plo(1):,plo(2):,plo(3):,plo(4):) !CAPITALIZATION INTENTIONAL!
  integer, intent(IN)  :: loFl(*)
  real,    intent(IN)  :: FLX(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: FLY(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: FLZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
  !! -----------------------------------------------------

  integer :: i,j,k
  real    :: invNewDens

  !! End of data declaration ***********************************************
#ifdef DEBUG_UHD
98 format(A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format(A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print *, "plo" ,plo(1:MDIM+1)
  print 98,"Uin" ,(plo(i),i=1,4)
  print 99,"Uin" ,(lbound(Uin ,i),ubound(Uin ,i),i=1,4)
  print 99,"Uout",(lbound(Uout,i),ubound(Uout,i),i=1,4)
  print*,'tileLim:',tileLimits
#endif

  if (hy_fluxCorrect) then
     call Driver_abortFlash("hy_hllUnsplit: flux correction is not implemented!")
  end if

  if (hy_useGravity) then
     call Driver_abortFlash("hy_hllUnsplit: support for gravity not implemented!")
  end if

  if (.NOT.hy_updateHydroFluxes) then
     return
  end if

  !! ************************************************************************
  !! Unsplit update for conservative variables from n to n+1 time step
  !  instead of  call hy_hllUnsplitUpdate(blockID,dt,dtOld,del,datasize,blkLimits, ...)
  !! this section starts the update
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           Uout(VELX_VAR,i,j,k) = Uin(VELX_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(VELY_VAR,i,j,k) = Uin(VELY_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(VELZ_VAR,i,j,k) = Uin(VELZ_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(ENER_VAR,i,j,k) = Uin(ENER_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(DENS_VAR,i,j,k) = Uin(DENS_VAR,i,j,k) + flX(HY_DENS_FLUX,i,j,k) - flX(HY_DENS_FLUX,i+1,j,k)
           if (NDIM > 1) Uout(DENS_VAR,i,j,k) = Uout(DENS_VAR,i,j,k) + flY(HY_DENS_FLUX,i,j,k) - flY(HY_DENS_FLUX,i,j+1,k)
           if (NDIM > 2) Uout(DENS_VAR,i,j,k) = Uout(DENS_VAR,i,j,k) + flZ(HY_DENS_FLUX,i,j,k) - flZ(HY_DENS_FLUX,i,j,k+1)

        end do
     end do
  end do
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + flX(HY_XMOM_FLUX,i,j,k) - flX(HY_XMOM_FLUX,i+1,j,k)
           Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + flX(HY_YMOM_FLUX,i,j,k) - flX(HY_YMOM_FLUX,i+1,j,k)
           Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + flX(HY_ZMOM_FLUX,i,j,k) - flX(HY_ZMOM_FLUX,i+1,j,k)
           Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + flX(HY_ENER_FLUX,i,j,k) - flX(HY_ENER_FLUX,i+1,j,k)
        end do
     end do
  end do
#if NDIM > 1
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + flY(HY_XMOM_FLUX,i,j,k) - flY(HY_XMOM_FLUX,i,j+1,k)
           Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + flY(HY_YMOM_FLUX,i,j,k) - flY(HY_YMOM_FLUX,i,j+1,k)
           Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + flY(HY_ZMOM_FLUX,i,j,k) - flY(HY_ZMOM_FLUX,i,j+1,k)
           Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + flY(HY_ENER_FLUX,i,j,k) - flY(HY_ENER_FLUX,i,j+1,k)
        end do
     end do
  end do
#endif
#if NDIM > 2
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + flZ(HY_XMOM_FLUX,i,j,k) - flZ(HY_XMOM_FLUX,i,j,k+1)
           Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + flZ(HY_YMOM_FLUX,i,j,k) - flZ(HY_YMOM_FLUX,i,j,k+1)
           Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + flZ(HY_ZMOM_FLUX,i,j,k) - flZ(HY_ZMOM_FLUX,i,j,k+1)
           Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + flZ(HY_ENER_FLUX,i,j,k) - flZ(HY_ENER_FLUX,i,j,k+1)
        end do
     end do
  end do
#endif
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           invNewDens = 1.0 / Uout(DENS_VAR,i,j,k)
           Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) * invNewDens
           Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) * invNewDens
           Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) * invNewDens
           Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) * invNewDens
        end do
     end do
  end do

  !! Correct energy if necessary
  !  instead of  call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

#ifdef EINT_VAR
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           invNewDens = 1.0 / Uout(DENS_VAR,i,j,k)
           Uout(EINT_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) &
                - 0.5 * dot_product(Uout(VELX_VAR:VELZ_VAR,i,j,k),Uout(VELX_VAR:VELZ_VAR,i,j,k))
        end do
     end do
  end do
#endif
     
End Subroutine hy_hllUpdateSolution
