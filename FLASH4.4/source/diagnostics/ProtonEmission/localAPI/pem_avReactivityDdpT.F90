!!****if* source/diagnostics/ProtonEmission/localAPI/pem_avReactivityDdpT
!!
!! NAME
!!
!!  pem_avReactivityDdpT
!!
!! SYNOPSIS
!!
!!  pem_avReactivityDdpT (real (in) :: TkeV)
!!
!! DESCRIPTION
!!
!!  Given the temperature in keV, this function evaluates the Maxwell-averaged thermal
!!  reactivity in cm^3/s for the fusion reaction D(d,p)T:
!!
!!                      D + D -> T (1.01 MeV) + p (3.03 MeV)
!!
!!  The fitting functional from Bosch and Hale is used (Eqs. 12,13,14 in Bosch, H.-S. and
!!  Hale, G.M., Improved formulas for fusion cross-sections and thermal reactivities, Nuclear
!!  Fusion, 32, 611-31 (1992)). The valid temperature range is between 0.2 and 100 keV.
!!  The Bosch-Hale fitting function for the Maxwell-averaged thermal reactivity <R> can be
!!  written in the following form (see: S. Atzeni and J. Meyer-Ter-Vehn, The Physics of
!!  Inertial Fusion, International Series of Monographs on Physics 125, Oxford Science
!!  Publications (2004), Chapter 1.4.3):
!!
!!              <R> = C1 * X^(-5/6) * Y^2 * exp (-3 * X^(1/3) * Y)
!!
!!  where:
!!                                     C2 * T
!!                     X = 1 -  ---------------------
!!                              1 + C3 * T + C5 * T^2
!!
!!                     Y = C0 * T^(-1/3)
!!
!!  and the fitting parameters are given by:
!!
!!      C0 | 6.2696             (more exact from Bosch paper: C0 = (Bg^2/4)^(1/3) = 6.26958467)
!!      C1 | 3.7212 x 10^(-16)  (more exact from Bosch paper: C1 = 3.721204039 x 10^(-16))
!!      C2 | 3.4127 x 10^(-3)
!!      C3 | 1.9917 x 10^(-3)
!!      C5 | 1.0506 x 10^(-5)
!!
!!  
!!  For computational efficiency and stability reasons, we form the quantity:
!!
!!                     Z = C0 * (X/T)^(1/3)
!!
!!  and reformulate the expression for <R> as:
!!
!!               <R> = C1 * Sqrt(1/X)^3 * Z^2 * exp (-3Z)
!!
!! ARGUMENTS
!!
!!  TkeV : the deuterium temperature in keV
!!
!! NOTES
!!        
!!  none
!!
!!***

real function pem_avReactivityDdpT (TkeV)

  implicit none

  real, intent (in) :: TkeV

  pem_avReactivityDdpT = 0.0

  return
end function pem_avReactivityDdpT
