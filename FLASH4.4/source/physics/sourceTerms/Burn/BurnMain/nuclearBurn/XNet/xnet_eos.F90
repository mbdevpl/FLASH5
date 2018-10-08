!***************************************************************************************************
! eos_helm.f90 10/18/17
! Interface to FLASH EoS.
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

Module xnet_eos
  Use Driver_interface, Only: Driver_abortFlash
  Use Eos_interface, Only: Eos
  Use Simulation_interface, Only: Simulation_mapStrToInt
  Use xnet_types, Only: dp
  Implicit None

#include "Eos.h"
#include "Flash.h"
#include "constants.h"

  Real(dp) :: eosData(EOS_NUM)
  Real(dp) :: massFrac(SPECIES_BEGIN:SPECIES_END)
  Logical :: eosMask(EOS_VARS+1:EOS_NUM)
  !$omp threadprivate(eosData,massFrac)

  Integer :: inuc2unk(NSPECIES) ! Index in UNK for each XNet species

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes the FLASH EoS interface.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: iheat, iscrn
    Use nuclear_data, Only: nname
    Implicit None
    Character(5) :: tmp_name
    Character(4) :: unk_name
    Integer :: inuc, iunk
    eosMask = .false.
    if ( iheat > 0 ) then
      eosMask(EOS_CV) = .true.
      eosMask(EOS_DET) = .true.
    end if
    if ( iscrn > 0 ) eosMask(EOS_ETA) = .true.
#ifdef EOS_DETAT
    if ( iheat > 0 .and. iscrn > 0 ) eosMask(EOS_DETAT) = .true.
#endif
    Do inuc = 1, NSPECIES
      tmp_name = adjustl(nname(inuc))
      unk_name = tmp_name(1:4)
      call Simulation_mapStrToInt(unk_name,iunk,MAPBLOCK_UNK)
      inuc2unk(inuc) = iunk
    EndDo

    Return
  End Subroutine eos_initialize

  Subroutine xnet_eos_interface(t9,rho,y,ye,cv,etae,detaedt9)
    Use nuclear_data, Only: aa
    Use xnet_constants, Only: avn, epmev
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag, lun_stdout
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:)

    ! Ouput variables
    Real(dp), Intent(out) :: ye, cv, etae, detaedt9

    ! Local variables
    Real(dp) :: ytot, abar, zbar, z2bar, zibar
    Integer :: ierr

    ! Calculate Ye
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    If ( iscrn > 0 .or. iheat > 0 ) Then

      ! Load input variables for the eos
      eosData(EOS_TEMP) = t9*1.0e9
      eosData(EOS_DENS) = rho
      massFrac(inuc2unk) = y*aa

      ! Call the eos
      call Eos(MODE_DENS_TEMP,1,eosData,massFrac=massFrac,mask=eosMask,diagFlag=ierr)

      ! Convert units from ergs/g to MeV/nucleon and K to GK
      etae = eosData(EOS_ETA)
      cv = eosData(EOS_CV)*1.0e9/epmev/avn
#ifdef EOS_DETAT
      detaedt9 = eosData(EOS_DETAT)*1.0e9
#else
      detaedt9 = 0.0
#endif
      if ( ierr > 0 ) then
        Write(lun_stdout,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,etae,detaedt9
        call Driver_abortFlash('[xnet_eos] Error: too many Newton-Raphson iterations in Eos')
      endif
    Else
      etae = 0.0
      detaedt9 = 0.0
      cv = 0.0
    EndIf
    If ( idiag >= 3 ) Write(lun_diag,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,etae,detaedt9

    Return
  End Subroutine xnet_eos_interface

  Subroutine eos_screen(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine calls the Helmholtz EOS with the input temperature, density and composition. It
    ! returns the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd, &
      & thbim2, twm2bi
    Use xnet_controls, Only: idiag, iheat, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:), etae, detaedt9

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, bkt, abar, zbar, z2bar, zibar
    Real(dp) :: sratio, ae, dsratiodeta

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
    Call salpeter_ratio(etae,sratio,dsratiodeta)
    ztilde = sqrt(z2bar + zbar*sratio)
    If ( iheat > 0 ) Then
      dztildedt9 = 0.5*zbar/ztilde * dsratiodeta*detaedt9
    Else
      dztildedt9 = 0.0
    EndIf

    ! Calculate plasma quantities
    bkt = bok*t9
    lambda0 = sqrt(4.0*pi*rho*avn*ytot) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*avn*rho*ye))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi)
    If ( idiag >= 3 ) Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen

  Subroutine salpeter_ratio(eta,ratio,dratiodeta)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-----------------------------------------------------------------------------------------------
    Use fd, Only: fdm1h, fd1h, fdm3h
    Use xnet_controls, Only: iheat
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: eta

    ! Output variables
    Real(dp), Intent(out) :: ratio, dratiodeta

    ! Local variables
    Real(dp) :: fermip, fermim
    Real(dp) :: dfmdeta, dfpdeta

    ! Calculate f_(-1/2) and f_(1/2)
    fermim = fdm1h(eta)
    fermip = fd1h(eta)

    ! Evalutate the Salpeter ratio
    ratio = 0.5 * fermim/fermip
    If ( iheat > 0 ) Then
      dfmdeta = -0.5 * fdm3h(eta)
      dfpdeta = +0.5 * fermim
      dratiodeta = ratio * (dfmdeta/fermim - dfpdeta/fermip)
    Else
      dratiodeta = 0.0
    EndIf

    Return
  End Subroutine salpeter_ratio

  Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: aa, zz, zz2, zzi
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(:)

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot

    ! Calculate abundance moments
    ytot  = sum(y(:))
    atot  = sum(y(:) * aa(:))
    ztot  = sum(y(:) * zz(:))
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = sum(y(:) * zz2(:)) / ytot
    zibar = sum(y(:) * zzi(:)) / ytot
    ye = ztot / atot
    If ( idiag >= 3 ) Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye

    Return
  End Subroutine y_moment

End Module xnet_eos
