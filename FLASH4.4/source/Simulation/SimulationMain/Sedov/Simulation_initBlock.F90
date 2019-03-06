!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: tileDesc  )
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  tileDesc -        describes the block to initialize
!!
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_minRhoInit     Density floor for initial condition
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(solnData,tileDesc)

  use Simulation_interface, ONLY: Simulation_computeAnalytical
  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
     &  sim_nProfile, sim_rProf, sim_vProf, sim_pProf, sim_pExp, sim_rhoProf, &
     &  sim_tInitial, sim_gamma, sim_expEnergy, sim_pAmbient, sim_rhoAmbient, &
     &  sim_useProfileFromFile, sim_profileInitial, &
     &  sim_smallX, sim_smallRho, sim_minRhoInit, sim_smallP, sim_rInit, &
     &  sim_smallT, &
     &  sim_nSubZones, sim_xCenter, sim_yCenter, sim_zCenter, sim_inSubzones, sim_inszd, &
     sim_threadBlockList, sim_threadWithinBlock
  use Grid_interface, ONLY : Grid_getCellCoords, &
                             Grid_getCellVolumes, &
                             Grid_subcellGeometry, &
                             Grid_getDeltas
  use Grid_tile, ONLY : Grid_tile_t 
  use ut_interpolationInterface
 
  implicit none

#include "constants.h"
#include "Flash.h"
  
  real,              pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc
  
  integer,parameter :: op = 2
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk, kat
  real     ::  drProf
  real,allocatable,dimension(:) :: rProf, vProf, rhoProf, pProf
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek, eint
  real     ::  dist
  real     ::  vSub, rhoSub, pSub, errIgnored

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real,allocatable :: cellVolumes(:, :, :)
  integer,dimension(LOW:HIGH,MDIM) :: tileLimits
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits
  integer,dimension(MDIM) :: axis

!!$  real     :: dvSub(0:sim_nSubZones-1,0:(sim_nSubZones-1)*K2D)
  real,allocatable :: dvSub(:,:)
  real     :: dvc, quotinv

  real :: deltas(1:MDIM)

  if (sim_useProfileFromFile) then
     ! lazy initialization - should already have been done from Simulation_init
     if (sim_tinitial > 0.0) call sim_scaleProfile(sim_tinitial)
  end if

!!$  if (.NOT. sim_useProfileFromFile .OR. sim_tinitial .LE. 0.0) then
  if (sim_tinitial .LE. 0.0) then
     allocate(rProf(sim_nProfile))
     allocate(vProf(sim_nProfile))
     allocate(rhoProf(sim_nProfile))
     allocate(pProf(sim_nProfile))

  !
  !  Construct the radial samples needed for the initialization.
  !
     diagonal = (sim_xMax-sim_xMin)**2
     diagonal = diagonal + K2D*(sim_yMax-sim_yMin)**2
     diagonal = diagonal + K3D*(sim_zMax-sim_zMin)**2
     diagonal = sqrt(diagonal)
  
     drProf = diagonal / (sim_nProfile-1)
  
     do i = 1, sim_nProfile
        rProf(i)   = (i-1) * drProf
     enddo
!!$  !
!!$  !  If t>0, use the analytic Sedov solution to initialize the
!!$  !  code.  Otherwise, just use a top-hat.
!!$  !
!!$
!!$     if (sim_tInitial .gt. 0.) then
!!$        call set_analytic_sedov (sim_nProfile, rProf, rhoProf, pProf, & 
!!$          vProf, sim_tInitial, sim_gamma, sim_expEnergy, & 
!!$          sim_pAmbient, sim_rhoAmbient)
!!$     else
        do i = 1, sim_nProfile
           rhoProf(i) = sim_rhoAmbient
           pProf(i)   = sim_pAmbient
           vProf(i)   = 0.
           if (rProf(i) .le. sim_rInit) pProf(i) = sim_pExp
        enddo
!!$     
!!$     endif
  end if                        !useProfileFromFile

  ! get the coordinate information for the current block
  tileLimits = tileDesc%limits
  grownTileLimits = tileDesc%grownLimits
    
  call Grid_getDeltas(tileDesc%level, deltas)

  ! Find a real difference between z's if problem is >= 3D
  if (NDIM > 2) then
     dzz = deltas(KAXIS)
  ! Otherwise this problem is <= 2D, so dzz is meaningless
  else
     dzz = 0.0
  endif

  ! Find a real difference between y's if problem is >= 2D
  if (NDIM > 1) then
     dyy = deltas(JAXIS)
  ! Otherwise this problem is <= 1D, so dyy is meaningless
  else
    dyy = 0.0
  endif

  dxx = deltas(IAXIS)

  if (sim_tinitial > 0.0) then

     call Simulation_computeAnalytical(solnData,tileDesc,sim_tinitial)

     !There is no parallel region in Grid_initDomain and so we use the
     !same thread within block code for both multithreading strategies.

     !$omp parallel if (sim_threadBlockList .or. sim_threadWithinBlock) &
     !$omp default(none) &
     !$omp shared(grownTileLimits,xCoord,yCoord,zCoord,tileDesc,&
     !$omp sim_inSubzones,sim_nSubZones,sim_rProf,sim_minRhoInit,sim_smallRho,sim_smallP,&
     !$omp sim_smallX,sim_pProf,sim_rhoProf,sim_vProf,sim_gamma,sim_inszd,&
     !$omp sim_smallT,&
     !$omp solnData, &
     !$omp sim_xCenter,sim_yCenter,sim_zCenter) &
     !$omp private(i,j,k,ii,jj,kk,n,sumRho,sumP,sumVX,sumVY,sumVZ,&
     !$omp xx,yy,zz,xDist,yDist,zDist,dist,distInv,jLo,jHi,frac,vel,axis,&
     !$omp rho,p,vx,vy,vz,ek,e,eint,kat)

#if NDIM == 3
     !$omp do schedule(static)
#endif
     do k = grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS)

#if NDIM == 2
        !$omp do schedule(static)
#endif
        do j = grownTileLimits(LOW, JAXIS), grownTileLimits(HIGH, JAXIS)

#if NDIM == 1
           !$omp do schedule(static)
#endif
           do i = grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH, IAXIS)
              if (NSPECIES > 0) then
                 solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
                 solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
              end if
              solnData(DENS_VAR,i,j,k)=max(solnData(DENA_VAR,i,j,k), sim_smallRho)
              solnData(PRES_VAR,i,j,k)=max(solnData(PRSA_VAR,i,j,k), sim_smallP)
              solnData(ENER_VAR,i,j,k)=    solnData(ENRA_VAR,i,j,k)
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k)=    solnData(EINA_VAR,i,j,k)
#endif
              solnData(GAME_VAR,i,j,k)=sim_gamma
              solnData(GAMC_VAR,i,j,k)=sim_gamma
              solnData(VELX_VAR,i,j,k)=    solnData(VLXA_VAR,i,j,k)
#ifdef VLYA_VAR
              solnData(VELY_VAR,i,j,k)=    solnData(VLYA_VAR,i,j,k)
#endif
#ifdef VLZA_VAR
              solnData(VELZ_VAR,i,j,k)=    solnData(VLZA_VAR,i,j,k)
#endif
              solnData(TEMP_VAR,i,j,k)=sim_smallT
#ifdef BDRY_VAR
              solnData(BDRY_VAR,i,j,k)=    -1.0
#endif
           enddo
#if NDIM == 1
           !$omp end do nowait
#endif
        enddo
#if NDIM == 2
        !$omp end do nowait
#endif
     enddo
#if NDIM == 3
     !$omp end do nowait
#endif
     !$omp end parallel


     RETURN                     ! DONE here!
  end if


  allocate(xCoord(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS))); xCoord = 0.0
  allocate(yCoord(grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS))); yCoord = 0.0
  allocate(zCoord(grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS))); zCoord = 0.0
  allocate(cellVolumes(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS), &
                       grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS), &
                       grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS)))

  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          xCoord)
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          yCoord)
#endif
#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          zCoord)
#endif
  call Grid_getCellVolumes(tileDesc%level, &
                           lbound(cellVolumes), ubound(cellVolumes), cellVolumes)

  !
  !     For each cell
  !  

  !There is no parallel region in Grid_initDomain and so we use the
  !same thread within block code for both multithreading strategies.

  !$omp parallel if (sim_threadBlockList .or. sim_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(grownTileLimits,xCoord,yCoord,zCoord,cellVolumes,tileDesc,&
  !$omp sim_inSubzones,sim_nSubZones,sim_rProf,sim_minRhoInit,sim_smallRho,sim_smallP,&
  !$omp sim_smallX,sim_pProf,sim_rhoProf,sim_vProf,sim_gamma,sim_inszd,&
  !$omp rProf,pProf,rhoProf,vProf,dxx,dyy,dzz,&
  !$omp sim_smallT,&
  !$omp sim_useProfileFromFile,sim_tinitial,errIgnored,solnData, &
  !$omp sim_rhoAmbient,sim_pAmbient, &
  !$omp sim_xCenter,sim_yCenter,sim_zCenter) &
  !$omp private(i,j,k,ii,jj,kk,n,sumRho,sumP,sumVX,sumVY,sumVZ,&
  !$omp xx,yy,zz,xDist,yDist,zDist,dist,distInv,jLo,jHi,frac,vel,axis,&
  !$omp rho,p,vx,vy,vz,ek,e,eint,kat,rhoSub,pSub,vSub,dvc,quotinv,dvSub)

  allocate(dvSub(0:sim_nSubZones-1,0:(sim_nSubZones-1)*K2D))

#if NDIM == 3
  !$omp do schedule(static)
#endif
  do k = grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS)
     zz = zCoord(k)

#if NDIM == 2
     !$omp do schedule(static)
#endif
     do j = grownTileLimits(LOW, JAXIS), grownTileLimits(HIGH, JAXIS)
        yy = yCoord(j)
        
#if NDIM == 1
        !$omp do schedule(static)
#endif
        do i = grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH, IAXIS)
           xx = xCoord(i)

           dvc = cellVolumes(i, j, k)
           call Grid_subcellGeometry(   sim_nSubZones, &
                                     1+(sim_nSubZones-1)*K2D, &
                                     1+(sim_nSubZones-1)*K3D, &
                                     dvc, dvSub, &
                                     xCoord(i)-0.5*dxx, xCoord(i)+0.5*dxx)

           sumRho = 0.
           sumP   = 0.
           sumVX  = 0.
           sumVY  = 0.
           sumVZ  = 0.
           
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + ((real(kk)+0.5)*sim_inSubzones-.5)*dzz 
              zDist = (zz - sim_zCenter) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + ((real(jj)+0.5)*sim_inSubzones-.5)*dyy
                 yDist = (yy - sim_yCenter) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + ((real(ii)+0.5)*sim_inSubzones-.5)*dxx
                    xDist = xx - sim_xCenter
                    
                    dist    = sqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
                    call sim_find (rProf, sim_nProfile, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We do a
                    !  linear interpolation of the quantities at jLo and jHi and sum those.
                    ! 
                    if (jLo .eq. 0) then
                       jLo = 1
                       jHi = 1
                       frac = 0.
                    else if (jLo .eq. sim_nProfile) then
                       jLo = sim_nProfile
                       jHi = sim_nProfile
                       frac = 0.
                    else
                       jHi = jLo + 1
                       frac = (dist - rProf(jLo)) / & 
                         (rProf(jHi)-rProf(jLo))
                    endif

                    pSub   =  pProf(jLo) + frac*(pProf(jHi)  - pProf(jLo))

                    rhoSub =  rhoProf(jLo) + frac*(rhoProf(jHi)- rhoProf(jLo))
                    rhoSub = max(rhoSub, sim_minRhoInit)

                    vSub   = vProf(jLo) + frac*(vProf(jHi)  - vProf(jLo))

                    ! 
                    !   Now total these quantities.   Note that  v is a radial velocity; 
                    !   we multiply by the tangents of the appropriate angles to get
                    !   the projections in the x, y and z directions.
                    !
                    sumP = sumP + pSub * dvSub(ii,jj)
                    
                    sumRho = sumRho + rhoSub * dvSub(ii,jj)
                    
                    vel = vSub * dvSub(ii,jj)
                    
                    sumVX  = sumVX  + vel*xDist*distInv
                    sumVY  = sumVY  + vel*yDist*distInv
                    sumVZ  = sumVZ  + vel*zDist*distInv
                    
                 enddo
              enddo
           enddo
           
!!$           quotinv = sim_inszd
           quotinv = 1.0 / dvc
           rho = max(sumRho * quotinv, sim_smallRho)
           p   = max(sumP   * quotinv, sim_smallP)
           vx  = sumVX  * quotinv
           vy  = sumVY  * quotinv
           vz  = sumVZ  * quotinv
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           e   = p/(sim_gamma-1.)
           eint= e/rho
           e   = e/rho + ek
           e   = max (e, sim_smallP)
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
              solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
           end if
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k)=eint
#endif
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
           solnData(TEMP_VAR,i,j,k)=sim_smallT
#ifdef BDRY_VAR
           solnData(BDRY_VAR,i,j,k)=    -1.0
#endif
        enddo
#if NDIM == 1
  !$omp end do nowait
#endif
     enddo
#if NDIM == 2
  !$omp end do nowait
#endif
  enddo
#if NDIM == 3
  !$omp end do nowait
#endif

  deallocate(dvSub)
  !$omp end parallel

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  deallocate(cellVolumes)

  deallocate(rProf)
  deallocate(vProf)
  deallocate(rhoProf)
  deallocate(pProf)

  return
end subroutine Simulation_initBlock




!  Subroutine:  set_analytic_sedov()

!  Description: Given a set of arrays to store radius, density, pressure, and
!               velocity, together with a time and a ratio of specific heats,
!               generate the analytical solution to the Sedov problem.

!               Currently this routine is limited to gamma=1.4 (I've hardwired
!               the dimensionless constant beta) and will complain if it gets
!               anything else.


subroutine set_analytic_sedov (N, r, rho, p, v, t, gamma, E, & 
     p_ambient, rho_ambient)

  !==============================================================================

  implicit none

  !  Arguments
  integer, intent(IN)     :: N
  real, intent(IN)        :: gamma, t, rho_ambient, E, p_ambient
  real, intent(IN), dimension(N)  :: r
  real, intent(OUT), dimension(N) :: rho, v, p

  !  Local variables

  real    beta, R0, dr, xi, nu1, nu2, nu3, nu4, nu5, VV, G, Z, & 
       kappa, zeta, epsilon, c2sqr, k, gamp1, gam7, gamm1
  integer i

  !==============================================================================

  if (gamma .ne. 1.4) then
     write (*,*) 'Warning!  Simulation_initBlock() found gamma<>1.4 and t>0.'
     write (*,*) '          Analytical initial conditions will be'
     write (*,*) '          wrong.  Assuming beta=1.033...'
  endif

  !               Compute dimensionless scaling constant and explosion radius.

  beta = 1.033
  R0   = beta * (E*t*t/rho_ambient)**0.2
  if (R0 == 0.0) then
     R0 = max(r(1),r(min(2,N))) * 1.e+30
  end if

  !               Compute exponents for self-similar solution.

  nu1 = - (13.*gamma*gamma - 7.*gamma + 12.) / & 
       ((3.*gamma - 1.) * (2.*gamma+1.))
  nu2 = 5. * (gamma - 1.) / (2.*gamma + 1.)
  nu3 = 3. / (2.*gamma + 1.)
  nu4 = - nu1 / (2. - gamma)
  nu5 = - 2. / (2. - gamma)

  !               Other useful combinations of gamma.

  gamp1 = gamma + 1.E0
  gamm1 = gamma - 1.E0
  gam7  = 7.E0 - gamma
  k     = gamp1 / gamm1

  !==============================================================================

  !               Generate the solution.

  do i = 1, N

     xi = r(i) / R0                ! Fraction of explosion radius.

     if (xi .le. 1.) then          ! Inside the explosion.

        ! Compute V(xi) using bisection.
        ! See Landau & Lifshitz for notation.
        call compute_sedov_v (xi, gamma, nu1, nu2, VV)

        G = k * (k*(gamma*VV-1.))**nu3 * & 
             (gamp1/gam7*(5.-(3.*gamma-1.)*VV))**nu4 * & 
             (k*(1.-VV))**nu5

        rho(i) = rho_ambient * G
        if (t.EQ.0.0) then
           v(i) = sqrt(HUGE(v(i))) * 1.e-10
           p(i) =      HUGE(p(i))  * 1.e-10
           CYCLE
        end if

        v(i)   = 2.*r(i)*VV / (5.*t)

        if (xi .le. 1.E-6) then     ! Use asymptotic r->0 solution.
           kappa = ( (0.5*gamp1/gamma)**2. * & 
                (gamp1/gam7* & 
                (5.-(3.*gamma-1.)/gamma))**nu1 )**(1./nu2) * & 
                gamm1/gamp1 / gamma

           epsilon = k**(nu5+1.) * (k*gamma*kappa)**nu3 * & 
                (gamp1/gam7*(3.*gamma-1.))**nu4 * & 
                ((2.*gamma+1)/gamma/(3.*gamma-1.)) * & 
                (gamm1/gamma)

           zeta = gamm1*gamm1/(2.*gamma*gamma*kappa)
           p(i) = rho_ambient/gamma * 0.16*(R0/t)**2 * epsilon*zeta
        else
           Z = gamma * gamm1 * (1.-VV) * VV**2 / (2.*(gamma*VV-1.))
           c2sqr = 0.16*(r(i)/t)**2 * Z
           p(i) = rho(i) * c2sqr / gamma
        endif

     else                          ! Outside the explosion.

        rho(i) = rho_ambient
        p(i)   = p_ambient
        v(i)   = 0.

     endif

  enddo

  return
end subroutine set_analytic_sedov



!  Subroutine:  compute_sedov_v()

!  Description: Compute the dimensionless velocity function V(xi) in the
!               Sedov problem using bisection.  See Landau & Lifshitz for
!               notation.

subroutine compute_sedov_v (xi, gamma, nu1, nu2, V)

  !==============================================================================

  implicit none

!  Arguments, LBR guessed intent on these
  real, intent(IN)  :: xi, gamma, nu1, nu2
  real, intent(OUT) :: V

! local variables
  real      VL, VR, xiL, xiR, Vmid, ximid, sedov_vfunc, tolerance, logxi
  integer   n_iter, n_iter_max
  parameter (n_iter_max = 500, tolerance = 1.E-6)

  !==============================================================================

  if (xi .le. 1.E-6) then         ! Use asymptotic xi->0 solution.

     !CD: The expression gamma*VV-1 on line 387 in set_analytic_sedov
     !(after this function call) should be 0.0 but generates a tiny value
     !(gdb) p gamma*VV-1.
     ! $1 = -1.1102230246251565e-16
     !This causes an FPE when it is raised to the power nu3.
     !When run in gdb it returns: 
     !"Cannot perform exponentiation: Numerical argument out of domain"
     !We can avoid the FPE by adding 1.E-8 to the original value of V:
     !V = 1./gamma
     V = 1./gamma + 1.E-8

  else                            ! Do bisection root search.


     logxi = alog(xi)             !was dlog but this fails on xlf compiler 
     VL = 1./gamma
     VR = 2./gamma
     xiL = sedov_vfunc(VL, gamma, nu1, nu2)
     xiR = sedov_vfunc(VR, gamma, nu1, nu2)
     n_iter = 1
10   Vmid = 0.5 * (VL + VR)
     ximid = sedov_vfunc(Vmid, gamma, nu1, nu2)
     if ((abs(ximid - logxi) .le. tolerance) .or. & 
          (n_iter .gt. n_iter_max)) goto 20
     n_iter = n_iter + 1
     if (ximid .gt. logxi) then
        VR = Vmid
     else
        VL = Vmid
     endif
     goto 10

#ifdef DEBUG_SIM
20   if (n_iter .gt. n_iter_max) & 
          write (*,*) 'compute_sedov_v:  did not reach ', & 
          'max precision for xi = ', xi
     V = Vmid
#else
20   V = Vmid
#endif

  endif

  return
end subroutine compute_sedov_v



!******************************************************************************

!  Function:    sedov_vfunc()

!  Description: Function to use in bisection root search (compute_sedov_v()).

real function sedov_vfunc (V, gamma, nu1, nu2)

  implicit none

! arguments, LBR guessed intent on these
  real, intent(IN)  :: V, gamma, nu1, nu2

!  local variables
  real gamp1, gamm1, gam7, k, xi, Vfpe


  !CD: The sub expression gamma*V-1 in xi expression causes a problem when 
  !V = 1./gamma.  Mathematically this is 0.0, but rounding gives a 
  !subexpression value of order -1.E-16.  This value is then passed
  !to the alog function (See alog(gamma*Vfpe-1.)) which gives -Infinity
  !for alog(0.0) and NaN for alog(-1.E-16).  Neither is good, so I 
  !guard against both cases by adding 1.E-8 to V.
  if ((gamma*V-1.) <= 0.0) then
     Vfpe = 1./gamma + 1.E-8
  else
     Vfpe = V
  end if
 

  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam7  = 7. - gamma
  k     = gamp1 / gamm1

  xi = nu1*alog(5.-(3.*gamma-1.)*Vfpe) + & 
       nu2*alog(gamma*Vfpe-1.) - & 
       nu1*alog(gam7/gamp1) - nu2*alog(gamm1/gamp1) - & 
       2.*alog(0.5*gamp1)
  sedov_vfunc = 0.2 * xi

  return
end function sedov_vfunc



!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, N, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: N
  integer, intent(OUT):: i
  real, intent(IN)    :: x(N), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(N)) then

     i = N

  else

     il = 1
     ir = N
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find
