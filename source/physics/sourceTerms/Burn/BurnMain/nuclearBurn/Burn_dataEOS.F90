!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_dataEOS
!!
!! NAME
!!  
!!  Burn_dataEOS
!!
!!
!! SYNOPSIS
!! 
!!  use Burn_dataEOS
!!
!! DESCRIPTION
!!
!!  Was eos_common.fh
!!
!! NOTES 
!!
!!  No idea where these are actually set, however!
!!
!!***

module Burn_dataEOS

!..equation of state communication:

!..btemp    = temperature
!..bden     = density
!..abar     = average number of nucleons per nuclei
!..zbar     = average number of protons per nuclei
!..z2bar    = square of zbar
!..ytot1    = total number of moles per gram
!..bye      = electron mole number


      real           ::  btemp,bden,abar,zbar,z2bar,ytot1,bye

!$omp threadprivate(btemp,bden,abar,zbar,z2bar,ytot1,bye)

end module Burn_dataEOS
