!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestRefine/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the Sod problem
!!
!!***

module Simulation_data
    implicit none

    integer, parameter :: MIN_REFINE_LEVEL = 1
    integer, parameter :: MAX_REFINE_LEVEL = 4
    integer, parameter :: NUM_LEVELS = MAX_REFINE_LEVEL - MIN_REFINE_LEVEL + 1 

    type blocks_t
      integer, allocatable :: blocks(:, :) 
    end type blocks_t
    
    ! Store the leaf blocks at each level
    type(blocks_t), save :: leaves(MIN_REFINE_LEVEL:MAX_REFINE_LEVEL)
end module Simulation_data

