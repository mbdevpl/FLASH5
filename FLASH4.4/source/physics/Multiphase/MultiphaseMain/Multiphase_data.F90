!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_data
!!
!! NAME
!!
!!  Multiphase_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Multiphase_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Multiphase module.
!!
!!***
 
 
module Multiphase_data

#include "Flash.h"
#include "constants.h"

  real, save :: mph_rho1
  real, save :: mph_rho2

  real, save :: mph_vis1
  real, save :: mph_vis2

  real, save :: mph_sten

  real, save :: mph_crmx, mph_crmn

  integer, save :: mph_lsit
  integer, save :: mph_inls

  integer, save :: mph_meshMe
  integer, save :: mph_meshNumProcs
  integer, save :: mph_meshComm


end module Multiphase_data
