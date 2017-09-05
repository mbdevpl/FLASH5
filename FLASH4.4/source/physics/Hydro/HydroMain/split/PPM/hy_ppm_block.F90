!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_block
!!
!! NAME
!!
!!  hy_ppm_block
!!
!! SYNOPSIS
!!
!!  call hy_ppm_block(integer(IN):: hy_meshMe,
!!                integer(IN):: blockID,
!!                integer(IN):: sweepDir, 
!!                real(IN)   :: dt, 
!!                real(IN)   :: dtOld, 
!!                integer(IN):: blkLimits(HIGH,MDIM),
!!                integer(IN):: blkLimitsGC(HIGH,MDIM),
!!                integer(IN):: bcs(2,MDIM),
!!                integer(IN):: numCells, 
!!                integer(IN):: numguard,
!!                real(IN)   :: primaryCoord(numCells),
!!                real(IN)   :: primaryLeftCoord(numCells),
!!                real(IN)   :: primaryRightCoord(numCells),
!!                real(IN)   :: primaryDx(numCells),
!!                real(IN)   :: secondCoord(numCells),
!!                real(IN)   :: thirdCoord(numCells),
!!                real(IN)   :: radialCoord(numCells),
!!                real(IN)   :: ugrid(numCells),
!!                real(OUT)  :: tempArea (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempGrav1d_o(blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempGrav1d  (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempDtDx    (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempFict    (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempAreaLeft(blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(OUT)  :: tempFlx(NFLUXES, blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real(IN)   :: shock       (blkLimitsGC(1,1):blkLimitsGC(2,1), blkLimitsGC(1,2):blkLimitsGC(2,2), blkLimitsGC(1,3):blkLimitsGC(2,3)),
!!                real, pointer :: solnData(:,:,:,:)
!!
!! DESCRIPTION
!!
!!  This routine takes a pointer to a block of grid data, picks out 1d
!!  slices, and applies the ppm algorithm, 1 row at a time, calculating
!!  fluxes.  The fluxes are stored and returned in the tempFlx/y/z
!!  arrays.  Other quantities used in the update step are calculated, as
!!  well, such as the geometry factors.
!!
!! ARGUMENTS
!!
!!   hy_meshMe -- my Processor Number
!!
!!   blockID -- My block number
!!
!!   sweepDir --  direction in which to do 1d hydro sweeps
!!
!!   dt --         timestep to advance through
!!
!!   dtOld --     previous timestep
!!
!!   blkLimits --  array holding upper and lower index limits of interior block cells (no GC)
!!
!!   blkLimitsGC --  array holding the upper and lower index limits of an entire block (including GC)
!!
!!   bcs --  boundary conditions for domain boundaries as obtained from Grid_getBlkBC
!!
!!   numCells --  tells this routine how long each 1d slice will be
!!
!!   numguard --  tells this routine how many guardcells surround the interior
!!
!!   primaryCoord --  positions of cell centers of the 1d slice
!!
!!   primaryLeftCoord --  positions of left interfaces
!!
!!   primaryRightCoord --  positions of right interfaces
!!
!!   primaryDx --  width of cells in the sweep direction
!!
!!   secondCoord --  for an x sweep: y coordinates; for a y sweep: x coord; for a z sweep: xcoord
!!
!!   thirdCoord --   for an x sweep: z coordinates; for a y sweep: z coord; for a z sweep: ycoord
!!
!!   radialCoord --  
!!
!!   ugrid --  this is for moving grid velocities
!!
!!   tempArea --  in case we're in non-cartesian geometry, compute cell face areas; used in the update step
!!
!!   tempGrav1d_o --  gravitational acceleration from the previous step, at cell edges
!!
!!   tempGrav1d --  gravitational acceleration from the current step 
!!
!!   tempDtDx --  dt/dx factor for updating quantities from fluxes
!!
!!   tempFict --  represents geometry related forces, e.g., centrifugal; used to update velocities 
!!
!!   tempFlx --   fluxes for all flux variables and all cells of the block; are used to update the solution, and
!!               given to Grid package for flux correction routines.
!!   tempAreaLeft --  cell face areas to the left
!!
!!
!!   shock --  0 if there is no shock or 1 if there is a shock for each cell; calculated by Hydro_detectShock 
!!
!!   solnData --  a pointer to the cell-centered data for the whole block
!!
!!***


! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx


#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif


subroutine hy_ppm_block( hy_meshMe,block,sweepDir, dt, dtOld, &
                         lim,limGC,bcs,  &
                         blkLimits,blkLimitsGC,  &
                         numCells,numguard, &
                         primaryCoord ,     &
                         primaryLeftCoord , &
                         primaryRghtCoord , &
                         primaryDx        , &
                         secondCoord      , &
                         thirdCoord       , &
                         radialCoord      , &
                         ugrid            , &
                         tempArea, tempGrav1d_o, tempGrav1d,          &
                         tempDtDx, tempFict, tempAreaLeft,            &
                         tempFlx,                & 
                         shock, solnData )


  use Hydro_data, ONLY : hy_useGravity, hy_dirGeom,&
                         hy_dela, &
                         hy_dp , hy_du , hy_dut, hy_dutt, &
                         hy_drho, hy_dgame, hy_dgamc, hy_dgrav, &
                         hy_p6, hy_u6 , hy_ut6 , hy_utt6 , &
                         hy_rho6 , hy_game6 , hy_gamc6 , hy_grav6 , &
                         hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, hy_pwcubic,&
                         hy_gravl, hy_gravr,   &
                         hy_clft  , hy_plft  , hy_uttlft,   &
                         hy_ulft  , hy_vlft  , hy_utlft , &
                         hy_crght , hy_prght , hy_vrght ,         &
                         hy_urght , hy_utrght, hy_uttrgt, &
                         hy_gmelft, hy_gmergt, &
                         hy_gmclft, hy_gmcrgt, hy_pstor,hy_nriem, &
                         hy_deint,hy_eint6,hy_eiLft, hy_eiRght,     &
                         hy_dxn,hy_xn6,hy_xnlft, hy_xnrght,     &
                         hy_updateHydroFluxes, &
                         hy_useCellAreasForFluxes, &
                         hy_gravMass

  use Grid_interface, ONLY : Grid_getSingleCellVol, Grid_getBlkData
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata, ONLY : block_metadata_t
  !! LOCAL - effective SCRATCH space
  !! shock and shock_multid are for shock detecting
  !! u -> utrt are for PPM.
  !! These are slices that we pass to PPM.  Optimized for vector
  !! machines years ago, PPM needs these slices in a very specific
  !! way.  This is costly.  

  implicit none
#include "constants.h"
#include "Flash.h"
#include "PPM.h"

  !! ------------
  !! ---- ARGUMENTS
  type(block_metadata_t),intent(IN) :: block
  integer, intent(IN) :: hy_meshMe
  integer, intent(IN) :: sweepDir
  real,    intent(IN) :: dt, dtOld
  integer, intent(IN) :: numCells,numguard
  integer, intent(IN),dimension(LOW:HIGH,MDIM) :: limGC,lim,bcs, blkLimits, blkLimitsGC
  
  real, intent(OUT), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),          &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::       &
       tempArea,       &
       tempGrav1d_o,   &
       tempGrav1d,     &
       tempDtDx,       &
       tempFict,       &
       tempAreaLeft
  real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
  
  real, intent(OUT), DIMENSION(NFLUXES,                    &
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
       tempFlx  
  real, intent(IN), DIMENSION(numCells) :: primaryCoord ,  &
       primaryLeftCoord , &
       primaryRghtCoord , &
       primaryDx        , &
       secondCoord      , &
       thirdCoord       , &
       radialCoord      , &
       ugrid
  real,    pointer :: solnData(:,:,:,:) 

  !! ------------

  integer :: i1, i2, j1, j2, k1, k2
  integer :: i, j, k, ii,ilo,ihi,jlo,jhi,klo,khi
  integer :: iloGc,ihiGc,jloGc,jhiGc,kloGc,khiGc
  real    :: xbot, xtop, ybot, ytop, ylft, yrgt, zlft, zrgt
  integer :: numIntCells
  integer,parameter :: numXn=NSPECIES+NMASS_SCALARS
  integer :: sp,istat
  ! for gravity accumulation, implemented by LBR 12/19/2006
  integer, DIMENSION(MDIM) :: point,len, startPos
  real               :: cellVolume
  real, allocatable :: faceAreas(:,:,:), cellVolumes(:,:,:)
  


  real, DIMENSION(numCells) :: dtdx, areaLeft, area, cvol, grav, ngrav, fict, shock_multid
  real, DIMENSION(numCells) :: rhoflx, uflx, pav, &
                               utflx, uttflx, eflx, eintflx,&
                               u, ut, utt, rho, p, e,  &
                               tmp, gamc, game,  &
                               uttp, utbt, utlt, utrt
  real, DIMENSION(numCells, NSPECIES+NMASS_SCALARS) :: xn, xnflx
  integer :: iglobal,jglobal,kglobal,iglobal1,iglobal2,jglobal1,jglobal2,kglobal1,kglobal2


  allocate(hy_dela(numCells),stat = istat) 
  if (istat==0)allocate(hy_dp(numCells),stat = istat) 
  if (istat==0)allocate( hy_du(numCells),stat = istat)  
  if (istat==0)allocate(hy_dut(numCells),stat = istat) 
  if (istat==0)allocate( hy_dutt(numCells),stat = istat) 
  if (istat==0)allocate( hy_drho(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgame(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgamc(numCells),stat = istat) 
  if (istat==0)allocate( hy_dgrav(numCells),stat = istat) 
  if (istat==0)allocate( hy_p6(numCells),stat = istat) 
  if (istat==0)allocate( hy_u6(numCells),stat = istat)
  if (istat==0)allocate( hy_ut6(numCells),stat = istat) 
  if (istat==0)allocate( hy_utt6(numCells),stat = istat) 
  if (istat==0)allocate( hy_rho6(numCells),stat = istat) 
  if (istat==0)allocate( hy_game6(numCells),stat = istat) 
  if (istat==0)allocate( hy_gamc6(numCells),stat = istat) 
  if (istat==0)allocate( hy_grav6(numCells),stat = istat) 
  if (istat==0)allocate( hy_pwl(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pwr(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_dpw(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pw6l(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pw6r(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_pwcubic(numCells),stat = istat) !! valid

  if (istat==0)allocate(hy_deint(numCells), stat = istat)
  if (istat==0)allocate(hy_eint6(numCells), stat = istat)
  if (istat==0)allocate(hy_eiLft(numCells),stat = istat) 
  if (istat==0)allocate( hy_eiRght(numCells),stat = istat) 

  if (istat==0)allocate(hy_dxn(numCells,numXN), stat = istat)
  if (istat==0)allocate(hy_xn6(numCells,numXN), stat = istat)
  if (istat==0)allocate(hy_xnlft(numCells,numXn),stat = istat) 
  if (istat==0)allocate( hy_xnrght(numCells,numXn),stat = istat) 

  if (istat==0)allocate(hy_gravl(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gravr(numCells),stat = istat) !! valid
  if (istat==0)allocate(hy_clft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_plft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_uttlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_ulft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_vlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_utlft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_crght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_prght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_vrght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_urght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_utrght(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_uttrgt(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmelft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmergt(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmclft(numCells),stat = istat) !! valid
  if (istat==0)allocate( hy_gmcrgt(numCells),stat = istat) !! valid
  

  if (istat==0)allocate( hy_pstor(hy_nriem + 2),stat = istat)
  if (istat .NE. 0) then
     call Driver_abortFlash("Memory allocation error in subroutine hy_block!")
  end if

!    Compute fluxes:  copy fluid data from 3D storage arrays        
!    into 1D work arrays, then operate on each row in turn.         
!    Update the solution in the 3D storage arrays, then save        
!    the computed fluxes for adjustment during the flux             
!    conservation step.                                             

!!$ IMPORTANT -- MOVE TO DATABASE. THIS IS A TEMPORARY HACK 
  grav(:) = 0.

  call Timers_start("hy_block")
! if we are using a hybrid Riemann solvers (i.e. using HLLE inside shocks), 
! then use a multi-dimensional shock detection -- this is more 
! accurate than the 1-d one done by the original PPM algorithm.
! Flags that signal in which cells a shock was detected are passed to this
! subroutine in the shock_multid array, if necessary.


  tempFlx(:,:,:,:) = 0.0
  select case (sweepDir)
  case (SWEEP_X)
     j1 = 1; j2 = 1
     k1 = 1; k2 = 1
     jglobal1 = 1; jglobal2 = 1
     kglobal1 = 1; kglobal2 = 1
     !! Loop over the interior to create the 1d slices
     iloGc = limGC(LOW,IAXIS)
     ihiGc = limGC(HIGH,IAXIS)
     ilo = lim(LOW,IAXIS)
     ihi = lim(HIGH,IAXIS)

     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)+1,     &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        startPos=1
        len(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
        len(IAXIS)=len(IAXIS)+1
        call Grid_getBlkData(block, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
             startPos, faceAreas, len)
        len(IAXIS)=len(IAXIS)-1
        allocate(cellVolumes(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        call Grid_getBlkData(block, CELL_VOLUME, 0, EXTERIOR, startPos, &
             cellVolumes,len)
     end if


     kglobal=lim(LOW,KAXIS)-1
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        point(KAXIS) = k
        kglobal=kglobal+1
        jglobal=lim(LOW,JAXIS)-1
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           point(JAXIS) = j
           jglobal=jglobal+1
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:blkLimits(HIGH,IAXIS)+1) = faceAreas(1:blkLimits(HIGH,IAXIS)+1,j,k)
              cvol(1:blkLimits(HIGH,IAXIS)) = cellVolumes(1:blkLimits(HIGH,IAXIS),j,k)
           end if
           shock_multid(1:len(IAXIS)) = shock(1:len(IAXIS),j,k)
           u(:)    = solnData( VELX_VAR,iloGc:ihiGc, jglobal,kglobal )
           ut(:)   = solnData( VELY_VAR, iloGc:ihiGc, jglobal,kglobal )
           utt(:)  = solnData( VELZ_VAR, iloGc:ihiGc, jglobal,kglobal )
           rho(:)  = solnData( DENS_VAR, iloGc:ihiGc, jglobal,kglobal )
           p(:)    = solnData( PRES_VAR, iloGc:ihiGc, jglobal,kglobal )
           e(:)    = solnData( ENER_VAR, iloGc:ihiGc, jglobal,kglobal )
           tmp(:)  = solnData( TEMP_VAR, iloGc:ihiGc, jglobal,kglobal )
           game(:) = solnData( GAME_VAR, iloGc:ihiGc, jglobal,kglobal )
           gamc(:) = solnData( GAMC_VAR, iloGc:ihiGc, jglobal,kglobal )

           if (numXN > 0) then
#ifdef INDEXREORDER
              !! Note: solnData is going to be flipped so this is OK even if it does not seem so
              xn(:,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iloGc:ihiGc,jglobal,kglobal)
#else
              xn(:,:) = transpose&
                (solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iloGc:ihiGc,jglobal,kglobal ))
#endif
           end if
           if (NDIM.ge.2) then
              j1 = j - 1 ; j2 = j + 1
              jglobal1=jglobal-1; jglobal2=jglobal+1
           end if
           if (NDIM.eq.3) then
              k1 = k - 1 ; k2 = k + 1
              kglobal1=kglobal-1; kglobal2=kglobal+1
           end if

           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(:) = solnData(VELY_VAR, iloGc:ihiGc, jglobal2, kglobal)
           utbt(:) = solnData(VELY_VAR, iloGc:ihiGc, jglobal1, kglobal)
           utrt(:) = solnData(VELZ_VAR, iloGc:ihiGc, jglobal, kglobal2)
           utlt(:) = solnData(VELZ_VAR, iloGc:ihiGc, jglobal, kglobal1)

           xbot = 0.
           xtop = 0.           
           ybot = secondCoord(j1)
           ytop = secondCoord(j2)
           ylft = 0.
           yrgt = 0.
           zlft = thirdCoord(k1)
           zrgt = thirdCoord(k2)
           numIntCells = ihi-ilo+1

!        call Timers_start("hydro_1d")
           call hydro_1d (block,numIntCells,numCells,numguard, bcs,&
                          sweepDir, hy_meshMe,dt, dtOld, &
                          j, k,                       &
                          hy_dirGeom(IAXIS), hy_useGravity,              &
                          xbot, xtop,                 &
                          ybot, ytop, ylft, yrgt,     &
                          zlft, zrgt, ugrid,          &
                          primaryCoord ,              &
                          primaryLeftCoord ,          &
                          primaryRghtCoord ,          &
                          primaryDx        ,          &
                          secondCoord      ,          &
                          thirdCoord       ,          &
                          radialCoord     ,           &
                          u, ut, utt, rho, p, e, tmp, game, gamc,   &
                          xn, utbt, uttp, utlt, utrt,               &
                          shock_multid,                             &
                          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                          rhoflx, uflx, pav, utflx, uttflx,         &
                          eflx, eintflx, xnflx)
!        call Timers_stop("hydro_1d")
           tempAreaLeft(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,j,k) = areaLeft(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1)
         if (hy_updateHydroFluxes) then
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              tempDtDx(i,j,k)     = dtdx(i)
              tempArea(i,j,k)     = area(i)
              tempGrav1d_o(i,j,k) = grav(i)
              tempGrav1d(i,j,k)   = ngrav(i)
              tempFict(i,j,k)     = fict(i)

              point(IAXIS) = i
              call Grid_getSingleCellVol(block,EXTERIOR,point,cellVolume)
              hy_gravMass(IAXIS) = hy_gravMass(IAXIS) + grav(i)*rho(i)*cellVolume
           enddo
           
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)+1
              tempFlx(RHO_FLUX,i,j,k)  = rhoflx(i)
              tempFlx(U_FLUX,i,j,k)    = uflx(i)
              tempFlx(P_FLUX,i,j,k)    = pav(i)
              tempFlx(UT_FLUX,i,j,k)   = utflx(i)
              tempFlx(UTT_FLUX,i,j,k)  = uttflx(i)
              tempFlx(E_FLUX,i,j,k)    = eflx(i)
              tempFlx(EINT_FLUX,i,j,k) = eintflx(i)
           enddo
           
           do sp = 1, numXn
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)+1
                 tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(i,sp)
              end do
           enddo
         else
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              tempDtDx(i,j,k)     = dtdx(i)
              tempArea(i,j,k)     = area(i)
              tempGrav1d_o(i,j,k) = grav(i)
              tempGrav1d(i,j,k)   = ngrav(i)
              tempFict(i,j,k)     = fict(i)
           end do
         end if
        
        end do !!j loop
     end do !!k loop

#if NDIM >= 2
  case (SWEEP_Y)
     k1 = 1; k2 = 1
     kglobal1 = 1; kglobal2 = 1
     !! Loop over the interior to create the 1d slices
     jlo = lim(LOW,JAXIS)
     jhi = lim(HIGH,JAXIS)
     jloGc = limGC(LOW,JAXIS)
     jhiGc = limGC(HIGH,JAXIS)
     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)+1,     &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        startPos=1
        len(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
        len(JAXIS)=len(JAXIS)+1
        call Grid_getBlkData(block, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
             startPos, faceAreas, len)
        len(JAXIS)=len(JAXIS)-1
        allocate(cellVolumes(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        call Grid_getBlkData(block, CELL_VOLUME, 0, EXTERIOR, startPos, &
             cellVolumes,len)
     end if
     kglobal=lim(LOW,KAXIS)-1
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        point(KAXIS) = k
        kglobal=kglobal+1
        iglobal=lim(LOW,IAXIS)-1
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           point(IAXIS) = i
           iglobal=iglobal+1
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:blkLimits(HIGH,JAXIS)+1) = faceAreas(i,1:blkLimits(HIGH,JAXIS)+1,k)
              cvol(1:blkLimits(HIGH,JAXIS)) = cellVolumes(i,1:blkLimits(HIGH,JAXIS),k)
           end if
           shock_multid(:) = shock(i,:,k)
           u(:)    = solnData( VELY_VAR,iglobal, jloGc:jhiGc, kglobal )
           ut(:)   = solnData( VELX_VAR,iglobal, jloGc:jhiGc, kglobal )
           utt(:)  = solnData( VELZ_VAR,iglobal, jloGc:jhiGc, kglobal )
           rho(:)  = solnData( DENS_VAR,iglobal, jloGc:jhiGc, kglobal )
           p(:)    = solnData( PRES_VAR,iglobal, jloGc:jhiGc, kglobal )
           e(:)    = solnData( ENER_VAR,iglobal, jloGc:jhiGc, kglobal )
           tmp(:)  = solnData( TEMP_VAR,iglobal, jloGc:jhiGc, kglobal )
           game(:) = solnData( GAME_VAR,iglobal, jloGc:jhiGc, kglobal )
           gamc(:) = solnData( GAMC_VAR,iglobal, jloGc:jhiGc, kglobal )
           if (numXN > 0) then
#ifdef INDEXREORDER
              xn(:,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iglobal,jloGc:jhiGc,kglobal)
#else
              xn(:,:) = transpose&
                &(solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1),iglobal,jloGc:jhiGc,kglobal ))
#endif
           end if
           i1 = i - 1 ; i2 = i + 1
           iglobal1 = iglobal - 1 ; iglobal2 = iglobal + 1
           if (NDIM.eq.3) then
              k1 = k - 1 ; k2 = k + 1
              kglobal1 = kglobal - 1 ; kglobal2 = kglobal + 1
           end if
           
           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(:) = solnData(VELX_VAR, iglobal2 ,jloGc:jhiGc, kglobal)
           utbt(:) = solnData(VELX_VAR, iglobal1,jloGc:jhiGc, kglobal)
           utrt(:) = solnData(VELZ_VAR,iglobal, jloGc:jhiGc, kglobal2)
           utlt(:) = solnData(VELZ_VAR,iglobal, jloGc:jhiGc, kglobal1)

           xbot = secondCoord(i1)
           xtop = secondCoord(i2)
           zlft = thirdCoord(k1)
           zrgt = thirdCoord(k2)
           numIntCells = jhi-jlo+1
!        call Timers_start("hydro_1d")
           call hydro_1d (block,numIntCells,numCells,numguard,bcs, &
                          sweepDir, hy_meshMe ,dt,dtOld, &
                          i, k,                       &
                          hy_dirGeom(JAXIS), hy_useGravity,                           &
                          xbot, xtop,                              &
                          ybot, ytop, ylft, yrgt,                  &
                          zlft, zrgt, ugrid,                       &
                          primaryCoord ,     &
                          primaryLeftCoord , &
                          primaryRghtCoord , &
                          primaryDx        , &
                          secondCoord      , &
                          thirdCoord       , &
                          radialCoord     , &
                          u, ut, utt, rho, p, e, tmp, game, gamc,  &
                          xn, utbt, uttp, utlt, utrt,              &
                          shock_multid,                            &
                          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                          rhoflx, uflx, pav, utflx, uttflx,        &
                          eflx, eintflx, xnflx)
!        call Timers_stop("hydro_1d")
           tempAreaLeft(i,blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,k) = areaLeft(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1)

         if (hy_updateHydroFluxes) then
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              tempDtDx(i,j,k)     = dtdx(j)
              tempArea(i,j,k)     = area(j)
              tempGrav1d_o(i,j,k) = grav(j)   !! Gradient in y direction
              tempGrav1d(i,j,k)   = ngrav(j)  !! Gradient in y direction + (dt/dtOld)*(difference in time)
              tempFict(i,j,k)     = fict(j)

              point(JAXIS) = j
              call Grid_getSingleCellVol(block,EXTERIOR,point,cellVolume)
              hy_gravMass(JAXIS) = hy_gravMass(JAXIS) + grav(j)*rho(j)*cellVolume !  grav(j)*rho(j)*cellVolume 

           enddo
!!$           
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              tempFlx(RHO_FLUX,i,j,k)  = rhoflx(j)
              tempFlx(U_FLUX,i,j,k)    = uflx(j)
              tempFlx(P_FLUX,i,j,k)    = pav(j)
              tempFlx(UT_FLUX,i,j,k)   = utflx(j)
              tempFlx(UTT_FLUX,i,j,k)  = uttflx(j)
              tempFlx(E_FLUX,i,j,k)    = eflx(j)
              tempFlx(EINT_FLUX,i,j,k) = eintflx(j)
           enddo
           
           do sp = 1, numXn
              do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
                 tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(j,sp)
              end do
           enddo
         else
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              tempDtDx(i,j,k)     = dtdx(j)
              tempArea(i,j,k)     = area(j)
              tempGrav1d_o(i,j,k) = grav(j)
              tempGrav1d(i,j,k)   = ngrav(j)
              tempFict(i,j,k)     = fict(j)
           end do
         end if
        end do !!j loop
     end do !!k loop

#endif
#if NDIM == 3 
  case (SWEEP_Z) 


     !! Loop over the interior to create the 1d slices
     klo = lim(LOW,KAXIS)
     khi = lim(HIGH,KAXIS)
     kloGc = limGC(LOW,KAXIS)
     khiGc = limGC(HIGH,KAXIS)

     if (hy_useCellAreasForFluxes) then
        allocate(faceAreas(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                           blkLimitsGC(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1))
        startPos=1
        len(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
        len(KAXIS)=len(KAXIS)+1
        call Grid_getBlkData(block, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
             startPos, faceAreas,len )
        len(KAXIS)=len(KAXIS)-1
        allocate(cellVolumes(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
        call Grid_getBlkData(block, CELL_VOLUME, 0, EXTERIOR, startPos,&
             cellVolumes, len)
     end if
     jglobal=lim(LOW,JAXIS)-1
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        point(JAXIS) = j
        jglobal=jglobal+1
        iglobal=lim(LOW,IAXIS)-1
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           point(IAXIS) = i
           iglobal=iglobal+1
           if (hy_useCellAreasForFluxes) then
              areaLeft(1:blkLimits(HIGH,KAXIS)+1) = faceAreas(i,j,1:blkLimits(HIGH,KAXIS)+1)
              cvol(1:blkLimits(HIGH,KAXIS)) = cellVolumes(i,j,1:blkLimits(HIGH,KAXIS))
           end if
           shock_multid(:) = shock(i,j,:)
           
           u(:)    = solnData( VELZ_VAR, iglobal, jglobal, kloGc:khiGc )
           ut(:)   = solnData( VELX_VAR, iglobal, jglobal, kloGc:khiGc )
           utt(:)  = solnData( VELY_VAR, iglobal, jglobal, kloGc:khiGc )
           rho(:)  = solnData( DENS_VAR, iglobal, jglobal, kloGc:khiGc )
           p(:)    = solnData( PRES_VAR, iglobal, jglobal, kloGc:khiGc )
           e(:)    = solnData( ENER_VAR, iglobal, jglobal, kloGc:khiGc )
           tmp(:)  = solnData( TEMP_VAR, iglobal, jglobal, kloGc:khiGc )
           game(:) = solnData( GAME_VAR, iglobal, jglobal, kloGc:khiGc )
           gamc(:) = solnData( GAMC_VAR, iglobal, jglobal, kloGc:khiGc )
           
           if (numXN > 0) then
#ifdef INDEXREORDER
              !! solnData will be flipped so this is OK
              xn(:,:) = solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1), iglobal, jglobal, kloGc:khiGc)
#else
              xn(:,:) = transpose&
                   &(solnData(SPECIES_BEGIN:(SPECIES_BEGIN+numXN-1), iglobal, jglobal, kloGc:khiGc))
#endif
           end if
           
           i1 = i - 1 ; i2 = i + 1
           j1 = j - 1 ; j2 = j + 1
           
           iglobal1 = iglobal - 1 ; iglobal2 = iglobal + 1
           jglobal1 = jglobal - 1 ; jglobal2 = jglobal + 1
           
           !!Initialise to zero (because we may not be overwriting the entire array)
           uttp = 0
           utbt = 0
           utrt = 0
           utlt = 0
           uttp(:) = solnData(VELX_VAR, iglobal2, jglobal, kloGc:khiGc)
           utbt(:) = solnData(VELX_VAR, iglobal1, jglobal, kloGc:khiGc)
           utrt(:) = solnData(VELY_VAR, iglobal, jglobal2, kloGc:khiGc)
           utlt(:) = solnData(VELY_VAR, iglobal, jglobal1, kloGc:khiGc)
           
           xbot = secondCoord(i1)
           xtop = secondCoord(i2)
           ylft = thirdCoord(j1)
           yrgt = thirdCoord(j2)
           numIntCells = khi-klo+1
           !        call Timers_start("hydro_1d")
           call hydro_1d (block,numIntCells,numCells,numguard,bcs, &
                sweepDir, hy_meshMe ,dt,dtOld, &
                i, j,                       &
                hy_dirGeom(KAXIS), hy_useGravity,        &
                xbot, xtop,                              &
                ybot, ytop, ylft, yrgt,                  &
                zlft, zrgt, ugrid,                       &
                primaryCoord ,     &
                primaryLeftCoord , &
                primaryRghtCoord , &
                primaryDx        , &
                secondCoord      , &
                thirdCoord       , &
                radialCoord     , &
                u, ut, utt, rho, p, e, tmp, game, gamc,  &
                xn, utbt, uttp, utlt, utrt,              &
                shock_multid,                            &
                dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
                rhoflx, uflx, pav, utflx, uttflx,        &
                eflx, eintflx, xnflx)
           !        call Timers_stop("hydro_1d")
           tempAreaLeft(i,j,blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) = areaLeft(blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
           
           if (hy_updateHydroFluxes) then
              do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
                 tempDtDx(i,j,k)     = dtdx(k)
                 tempArea(i,j,k)     = area(k)
                 tempGrav1d_o(i,j,k) = grav(k)
                 tempGrav1d(i,j,k)   = ngrav(k)
                 tempFict(i,j,k)     = fict(k)
                 
                 point(KAXIS) = k
                 call Grid_getSingleCellVol(block,EXTERIOR,point,cellVolume)
                 hy_gravMass(KAXIS) = hy_gravMass(KAXIS) + grav(k)*rho(k)*cellVolume
                 
              enddo
              
              do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1 
                 tempFlx(RHO_FLUX,i,j,k)  = rhoflx(k)
                 tempFlx(U_FLUX,i,j,k)    = uflx(k)
                 tempFlx(P_FLUX,i,j,k)    = pav(k)
                 tempFlx(UT_FLUX,i,j,k)  = utflx(k)
                 tempFlx(UTT_FLUX,i,j,k)  = uttflx(k)
                 tempFlx(E_FLUX,i,j,k)    = eflx(k)
                 tempFlx(EINT_FLUX,i,j,k) = eintflx(k)
              enddo
              
              do sp = 1, numXn
                 do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1  
                    tempFlx(SPECIES_FLUX_BEGIN+sp-1,i,j,k) = xnflx(k,sp)
                 end do
              enddo
              
           else
              do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
                 tempDtDx(i,j,k)     = dtdx(k)
                 tempArea(i,j,k)     = area(k)
                 tempGrav1d_o(i,j,k) = grav(k)
                 tempGrav1d(i,j,k)   = ngrav(k)
                 tempFict(i,j,k)     = fict(k)
              end do
           end if
           
        end do
     end do
     
#endif
  
  end select
  
  if (hy_useCellAreasForFluxes) then
     deallocate(cellVolumes)
     deallocate(faceAreas)
  end if
  deallocate( hy_pstor)
  
  deallocate(hy_dela) 
  deallocate(hy_dp) 
  deallocate( hy_du)  
  deallocate(hy_dut) 
  deallocate( hy_dutt) 
  deallocate( hy_drho) 
  deallocate( hy_dgame) 
  deallocate( hy_dgamc) 
  deallocate( hy_dgrav) 
  deallocate( hy_p6) 
  deallocate( hy_u6)
  deallocate( hy_ut6) 
  deallocate( hy_utt6) 
  deallocate( hy_rho6) 
  deallocate( hy_game6) 
  deallocate( hy_gamc6) 
  deallocate( hy_grav6) 
  deallocate( hy_pwl) 
  deallocate( hy_pwr) 
  deallocate( hy_dpw) 
  deallocate( hy_pw6l) 
  deallocate( hy_pw6r) 
  deallocate( hy_pwcubic) 
  
  deallocate(hy_deint)
  deallocate(hy_eint6)
  deallocate(hy_eiLft)
  deallocate( hy_eiRght)
  
  deallocate(hy_dxn)
  deallocate(hy_xn6)
  deallocate(hy_xnlft)
  deallocate( hy_xnrght)
  
  deallocate(hy_gravl) 
  deallocate( hy_gravr) 
  deallocate(hy_clft) 
  deallocate( hy_plft) 
  deallocate( hy_uttlft) 
  deallocate( hy_ulft) 
  deallocate( hy_vlft) 
  deallocate( hy_utlft) 
  deallocate( hy_crght) 
  deallocate( hy_prght) 
  deallocate( hy_vrght) 
  deallocate( hy_urght) 
  deallocate( hy_utrght) 
  deallocate( hy_uttrgt) 
  deallocate( hy_gmelft) 
  deallocate( hy_gmergt) 
  deallocate( hy_gmclft) 
  deallocate( hy_gmcrgt) 
  
  
  call Timers_stop("hy_block")
end subroutine hy_ppm_block

