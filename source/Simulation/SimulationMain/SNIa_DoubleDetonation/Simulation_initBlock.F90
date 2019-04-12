!!
!! Adapted from Dean M. Townsley (2009) SNIa_ddt setup
!!
!! This initializes the white dwarf and flame for Type Ia simulation
!!

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data
  use sim_local_interface, ONLY : sim_interpolate1dWd
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getCellCoords, Grid_getGeometry, &
     Grid_renormAbundance
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "constants.h"
#include "Flash.h"

  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc

  logical, parameter :: useGuardCell = .true.
  integer, dimension(2,MDIM) :: tileLimits, tileLimitsGC
  integer :: sizeX,sizeY,sizeZ,lev
  real, allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real, allocatable, dimension(:) :: xLeft, yLeft, zLeft
  real, allocatable, dimension(:) :: xRight, yRight, zRight
  real, dimension(MDIM) :: delta
  integer :: meshGeom

  real, allocatable, dimension(:) :: xinitial
  real :: radCenter, thtCenter, radCenterVol, dVol, ign_dist, dx, dy, dz, dr
  real :: radMin, radMax, x2Min, x2Max, y2Min, y2Max, z2Min, z2Max
  real :: velx, vely, velz, dens, temp, pres, eint, xsum, xmissing
  integer :: i, j, k, n, nunk

!==============================================================================

  call Grid_getGeometry(meshGeom)

  allocate(xinitial(sim_wd_nspec))

  ! Get the indices of the tileDescs
  tileLimits = tileDesc%limits
  tileLimitsGC = tileDesc%grownLimits
  lev = tileDesc%level

  allocate(xCenter(tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)))
  allocate(xRight(tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)))
  allocate(xLeft(tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)))
  allocate(yCenter(tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS)))
  allocate(yRight(tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS)))
  allocate(yLeft(tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS)))
  allocate(zCenter(tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)))
  allocate(zRight(tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)))
  allocate(zLeft(tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)))

  sizeX=SIZE(xCenter)
  sizeY=SIZE(yCenter)
  sizeZ=SIZE(zCenter)

  call Grid_getDeltas(tileDesc%level, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

  if (NDIM > 2) then
     call Grid_getCellCoords (KAXIS, CENTER, lev, tileLimits(LOW,:), tileLimits(HIGH,:), zCenter)
     call Grid_getCellCoords (KAXIS, LEFT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), zLeft)
     call Grid_getCellCoords (KAXIS, RIGHT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), zRight)
  end if
  if(NDIM>1) then
     call Grid_getCellCoords (JAXIS, CENTER, lev, tileLimits(LOW,:), tileLimits(HIGH,:), yCenter)
     call Grid_getCellCoords (JAXIS, LEFT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), yLeft)
     call Grid_getCellCoords (JAXIS, RIGHT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), yRight)
  end if
  call Grid_getCellCoords (IAXIS, CENTER, lev, tileLimits(LOW,:), tileLimits(HIGH,:), xCenter)
  call Grid_getCellCoords (IAXIS, LEFT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), xLeft)
  call Grid_getCellCoords (IAXIS, RIGHT_EDGE, lev, tileLimits(LOW,:), tileLimits(HIGH,:), xRight)

  do k = tileLimits(LOW,KAXIS), tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS), tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS), tileLimits(HIGH,IAXIS)

           !-----------------------------------------------
           !  determine state of material at this radius if unburned from external 1-d hyrdostatic model
           !-----------------------------------------------

           if ( meshGeom == SPHERICAL ) then
              radCenter = xCenter(i)
              thtCenter = yCenter(j)
              dVol = (4.0*PI/3.0)*(xRight(i)-xLeft(i))*(3.0*xLeft(i)*xRight(i) + (xRight(i)-xLeft(i))**2)
              radCenterVol = (4.0*PI/3.0)*xLeft(i)**3 + 0.5*dVol

              radMin = xLeft(i)
              radMax = xRight(i)
           else if ( meshGeom == CYLINDRICAL ) then
              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2)
              if ( yCenter(j) /= 0.0 ) then
                 thtCenter = atan( xCenter(i) / yCenter(j) )
              else
                 thtCenter = 0.5 * PI
              end if
              if ( thtCenter < 0.0 ) then
                 thtCenter = thtCenter + PI
              end if
              radCenterVol = (4.0*PI/3.0) * radCenter**3

              x2Min = xLeft(i)**2
              x2Max = xRight(i)**2

              if ( yLeft(j)*yRight(j) > 0.0 ) then
                 y2Min = min(yLeft(j)**2,yRight(j)**2)
              else
                 y2Min = 0.0
              end if
              y2Max = max(yLeft(j)**2,yRight(j)**2)

              radMin = sqrt( x2Min + y2Min )
              radMax = sqrt( x2Max + y2Max )
           else if ( meshGeom == CARTESIAN ) then
              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2 + zCenter(k)**2)
              thtCenter = acos( zCenter(k) / radCenter )
              radCenterVol = (4.0*PI/3.0) * radCenter**3

              if ( xLeft(i)*xRight(i) > 0.0 ) then
                 x2Min = min(xLeft(i)**2,xRight(i)**2)
              else
                 x2Min = 0.0
              end if
              x2Max = max(xLeft(i)**2,xRight(i)**2)

              if ( yLeft(j)*yRight(j) > 0.0 ) then
                 y2Min = min(yLeft(j)**2,yRight(j)**2)
              else
                 y2Min = 0.0
              end if
              y2Max = max(yLeft(j)**2,yRight(j)**2)

              if ( zLeft(k)*zRight(k) > 0.0 ) then
                 z2Min = min(zLeft(k)**2,zRight(k)**2)
              else
                 z2Min = 0.0
              end if
              z2Max = max(zLeft(k)**2,zRight(k)**2)

              radMin = sqrt( x2Min + y2Min + z2Min )
              radMax = sqrt( x2Max + y2Max + z2Max )
           else
              call Driver_abortFlash("Geometry not supported")
           end if

           call sim_interpolate1dWd(radCenterVol, radMin, radMax, dens, temp, xinitial)

           if(sim_useShell) then
             ! add a shell/belt
              if ( radCenter >= sim_radShellMin .and. radCenter <= sim_radShellMax .and. &
                 & thtCenter >= sim_thtShellMin .and. thtCenter <= sim_thtShellMax ) then
                 xinitial(:) = 0.0
                 if ( sim_wd_unk2spec(HE4_SPEC) > 0 ) xinitial(sim_wd_unk2spec(HE4_SPEC)) = sim_xhe4Shell
                 if ( sim_wd_unk2spec(C12_SPEC) > 0 ) xinitial(sim_wd_unk2spec(C12_SPEC)) = sim_xc12Shell
                 if ( sim_wd_unk2spec(NI56_SPEC) > 0 ) xinitial(sim_wd_unk2spec(NI56_SPEC)) = sim_xni56Shell
                 if ( sim_wd_unk2spec(O16_SPEC) > 0 ) xinitial(sim_wd_unk2spec(O16_SPEC)) = 1.0-sim_xhe4Shell-sim_xc12Shell-sim_xni56Shell

                 if ( dens > sim_densFluff ) then
                    dens = sim_densShellMult*dens
                    temp = sim_tempShellMult*temp
                 else
                    if ( sim_densShell > 0.0 ) then
                       dens = sim_densShell
                    end if
                    if ( sim_tempShell > 0.0 ) then
                       temp = sim_tempShell
                    end if
                 end if

              end if
           end if

           if (sim_ignite) then
              !-----------------------------------------------
              ! initialize burned region
              !-----------------------------------------------
              ! default to a spherical region centered at specified coordinates
              ! distance from center of ignition region
              ign_dist = (xCenter(i) - sim_ignX)**2
              if ( NDIM >= 2 ) ign_dist = ign_dist + (yCenter(j) - sim_ignY)**2
              if ( NDIM == 3 ) ign_dist = ign_dist + (zCenter(k) - sim_ignZ)**2
              ign_dist = sqrt(ign_dist)

              ! heat to ignition temp if zone center is within half zone-width of match
              if ( ign_dist <= sim_ignROuter + 0.5*dx .and. ign_dist >= sim_ignRInner - 0.5*dx ) then
                 temp = sim_ignTOuter
              end if

           end if ! sim_ignite

           ! Giant traffic cone
           if ( dens < sim_smallrho .or. temp < sim_smallt) then
              write(*,'(a,i5,a,3i3)') '[Before EOS] Bad value(s) on PE=',sim_globalMe,', (i,j,k)=',i,j,k
              write(*,'(a,2es15.7)') '  radCenter, thtCenter = ', radCenter, thtCenter
              write(*,'(a,2es15.7)') '     radMin,    radMax = ', radMin, radMax
              write(*,'(a,1es15.7)') '          radCenterVol = ', radCenterVol
              write(*,'(a,3es15.7)') '        x,y,z (center) = ', xCenter(i), yCenter(j), zCenter(k)
              write(*,'(a,3es15.7)') '        x,y,z   (left) = ', xLeft(i),   yLeft(j),   zLeft(k)
              write(*,'(a,3es15.7)') '        x,y,z  (right) = ', xRight(i),  yRight(j),  zRight(k)
              write(*,'(a,3es15.7)') '              dx,dy,dz = ', dx, dy, dz
              write(*,'(a,2es15.7)') '             dens,temp = ', dens, temp
              call Driver_abortFlash("[Simulation_initBlock] Bad values BEFORE EOS call.")
           end if

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           velx = 0.0
           vely = 0.0
           velz = 0.0

           solnData(VELX_VAR,i,j,k) = velx
           solnData(VELY_VAR,i,j,k) = vely
           solnData(VELZ_VAR,i,j,k) = velz
           solnData(DENS_VAR,i,j,k) = dens
           solnData(TEMP_VAR,i,j,k) = temp

           ! Copy interpolated values into UNK, do the renormalization later in Grid_renormAbundance
           solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = sim_smallx
           xsum = 0.0
           do n = 1, sim_wd_nspec
              nunk = sim_wd_spec2unk(n)
              if ( nunk /= NONEXISTENT ) then
                 solnData(nunk,i,j,k) = max(sim_smallx, min(1.0,xinitial(n)) )
                 xsum = xsum + solnData(nunk,i,j,k)
              end if
           end do

        end do
     end do
  end do

  call Grid_renormAbundance(tileDesc,tileLimits,solnData)

  call Eos_wrapped(MODE_DENS_TEMP,tileLimits,solnData)

  ! Giant traffic cone
  do k = tileLimits(LOW,KAXIS), tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS), tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS), tileLimits(HIGH,IAXIS)
           dens = solnData(DENS_VAR,i,j,k)
           temp = solnData(TEMP_VAR,i,j,k)
           pres = solnData(PRES_VAR,i,j,k)
           eint = solnData(EINT_VAR,i,j,k)
           if ( dens < sim_smallrho .or. temp < sim_smallt .or. pres < sim_smallp .or. eint < sim_smalle ) then
              write(*,'(a,i5,a,3i3)') '[After EOS] Bad value(s) on PE=',sim_globalMe,', (i,j,k)=',i,j,k
              write(*,'(a,3es15.7)') '        x,y,z (center) = ', xCenter(i), yCenter(j), zCenter(k)
              write(*,'(a,3es15.7)') '        x,y,z   (left) = ', xLeft(i),   yLeft(j),   zLeft(k)
              write(*,'(a,3es15.7)') '        x,y,z  (right) = ', xRight(i),  yRight(j),  zRight(k)
              write(*,'(a,3es15.7)') '              dx,dy,dz = ', dx, dy, dz
              write(*,'(a,4es15.7)') '   dens,temp,pres,eint = ', dens, temp, pres, eint
              call Driver_abortFlash("[Simulation_initBlock] Bad values AFTER EOS call.")
           end if
        end do
     end do
  end do

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(yCenter)
  deallocate(zLeft)
  deallocate(zRight)
  deallocate(zCenter)

  deallocate(xinitial)
  return
end subroutine Simulation_initBlock
