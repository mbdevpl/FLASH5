!!****if* source/diagnostics/ProtonEmission/localAPI/pem_avReactivityHe3dpHe4
!!
!! NAME
!!
!!  pem_avReactivityHe3dpHe4
!!
!! SYNOPSIS
!!
!!  pem_avReactivityHe3dpHe4 (real (in) :: TkeV)
!!
!! DESCRIPTION
!!
!!  Given the temperature in keV, this function evaluates the Maxwell-averaged thermal
!!  reactivity in cm^3/s for the fusion reaction He3(d,p)He4:
!!
!!                      He3 + D -> He4 (3.6 MeV) + p (14.7 MeV)
!!
!!  The fitting functional from Bosch and Hale is used (Eqs. 12,13,14 in Bosch, H.-S. and
!!  Hale, G.M., Improved formulas for fusion cross-sections and thermal reactivities, Nuclear
!!  Fusion, 32, 611-31 (1992)). The valid temperature range is between 0.5 and 190 keV.
!!  The Bosch-Hale fitting function for the Maxwell-averaged thermal reactivity <R> can be
!!  written in the following form (see: S. Atzeni and J. Meyer-Ter-Vehn, The Physics of
!!  Inertial Fusion, International Series of Monographs on Physics 125, Oxford Science
!!  Publications (2004), Chapter 1.4.3):
!!
!!              <R> = C1 * X^(-5/6) * Y^2 * exp (-3 * X^(1/3) * Y)
!!
!!  where:
!!                                C2 * T + C4 * T^2
!!                     X = 1 -  ---------------------
!!                              1 + C3 * T + C5 * T^2
!!
!!                     Y = C0 * T^(-1/3)
!!
!!  and the fitting parameters are given by:
!!
!!      C0 | 10.572
!!      C1 | 1.5116 x 10^(-14)
!!      C2 | 6.4192 x 10^(-3)
!!      C3 |-2.0290 x 10^(-3)
!!      C4 |-1.9108 x 10^(-5)
!!      C5 | 1.3578 x 10^(-4)
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

real function pem_avReactivityHe3dpHe4 (TkeV)

  implicit none

  real, intent (in) :: TkeV

  pem_avReactivityHe3dpHe4 = 0.0

  return
end function pem_avReactivityHe3dpHe4
