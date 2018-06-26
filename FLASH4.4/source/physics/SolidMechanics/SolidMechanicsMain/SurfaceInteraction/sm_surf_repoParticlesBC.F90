!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_repoParticlesBC
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_surf_repoParticlesBC(ibd)

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_sbData, ONLY : gr_sbBodyInfo
  use SolidMechanics_data, ONLY : sm_BodyInfo
  use Grid_data, ONLY  : gr_globalDomain, gr_domainBC
  implicit none

  ! IO variables
  integer, intent(in) :: ibd

  integer :: i,totpart
  real :: Lx,Ly,Lz

  integer ::nwx,nwy,nwz

  Lx = gr_globalDomain(HIGH,IAXIS) - gr_globalDomain(LOW,IAXIS)
  Ly = gr_globalDomain(HIGH,JAXIS) - gr_globalDomain(LOW,JAXIS)
#if NDIM == MDIM
  Lz = gr_globalDomain(HIGH,KAXIS) - gr_globalDomain(LOW,KAXIS)
#endif

  ! Initialize test array for body surface crossing boundary:
  sm_BodyInfo(ibd)%OnBoundary = .FALSE.

  totpart = gr_sbBodyInfo(ibd)%totalPart 
  do i=1,totpart  

     ! X Direction:
     if (gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) .lt. gr_globalDomain(LOW,IAXIS)) then

        ! This Body is OnBoundary at Low - x:
        sm_BodyInfo(ibd)%OnBoundary(LOW,IAXIS) = .TRUE.

        if (gr_domainBC(LOW,IAXIS)==PERIODIC) then

           nwx = int((-gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)+gr_globalDomain(LOW,IAXIS))/Lx)

           ! Add Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)+real(nwx+1)*Lx
        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with X non-periodic BC")   
        endif

     elseif (gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) .gt. gr_globalDomain(HIGH,IAXIS)) then

        ! This Body is OnBoundary at High - x:
        sm_BodyInfo(ibd)%OnBoundary(HIGH,IAXIS) = .TRUE.

        if (gr_domainBC(HIGH,IAXIS)==PERIODIC) then

           nwx = int((gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)-gr_globalDomain(HIGH,IAXIS))/Lx)

           ! Substract Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)-real(nwx+1)*Lx
        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with X non-periodic BC")   
        endif
     endif

     ! Y Direction:
     if (gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i) .lt. gr_globalDomain(LOW,JAXIS)) then

        ! This Body is OnBoundary at Low - y:
        sm_BodyInfo(ibd)%OnBoundary(LOW,JAXIS) = .TRUE.

        if (gr_domainBC(LOW,JAXIS)==PERIODIC) then

           nwy = int((-gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)+gr_globalDomain(LOW,JAXIS))/Ly)

           ! Add Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)+real(nwy+1)*Ly
        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with Y non-periodic BC")   
        endif

     elseif (gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i) .gt. gr_globalDomain(HIGH,JAXIS)) then

        ! This Body is OnBoundary at High - y:
        sm_BodyInfo(ibd)%OnBoundary(HIGH,JAXIS) = .TRUE.

        if (gr_domainBC(HIGH,JAXIS)==PERIODIC) then

           nwy = int((gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)-gr_globalDomain(HIGH,JAXIS))/Ly)

           ! Substract Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)-real(nwy+1)*Ly
        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with Y non-periodic BC")   
        endif
     endif

#if NDIM == MDIM
     ! Z Direction:
     if (gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i) .lt. gr_globalDomain(LOW,KAXIS)) then

        ! This Body is OnBoundary at Low - z:
        sm_BodyInfo(ibd)%OnBoundary(LOW,KAXIS) = .TRUE.

        if (gr_domainBC(LOW,KAXIS)==PERIODIC) then

           nwz = int((-gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)+gr_globalDomain(LOW,KAXIS))/Lz)

           ! Add Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)+real(nwz+1)*Lz

        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with Z non-periodic BC")   
        endif

     elseif (gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i) .gt. gr_globalDomain(HIGH,KAXIS)) then

        ! This Body is OnBoundary at High - z:
        sm_BodyInfo(ibd)%OnBoundary(HIGH,KAXIS) = .TRUE.

        if (gr_domainBC(HIGH,KAXIS)==PERIODIC) then

           nwz = int((gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)-gr_globalDomain(HIGH,KAXIS))/Lz)

           ! Substract Periodicity length:
           gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i) = gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)-real(nwz+1)*Lz
        else
           call Driver_abortFlash("sm_surf_repoParticlesBC: Particle outside domain with Z non-periodic BC")   
        endif
     endif
#endif

  enddo

  return

end subroutine sm_surf_repoParticlesBC
