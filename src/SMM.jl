
"""
# SMM.jl
This is the documentation for SMM.jl

"""

module SMM

  # Dependencies
  #-------------
  using BlackBoxOptim
  using JLD2
  using Plots
  using CSV
  using ProgressMeter
  using Calculus
  using DataStructures

  # Exports from BlackBoxOptim
  #---------------------------
  export best_candidate

  # Includes
  #---------
  # Types
  #------
  include("types.jl");

  # API
  #----
  include("api.jl")

  # General functions, useful at several places
  #---------------------------------------------
  include("generic.jl")

  # Functions to load and save
  #--------------------------
  include("save_load.jl")

  # Functions to minimize the objective function
  #---------------------------------------------
  include("optimize.jl")

  # Exports
  #--------
  # Functions and types in types.jl
  #----------------------------------
  export SMMOptions, SMMProblem
  export default_function, rosenbrock2d

  # Functions and types in api.jl
  #-------------------------------

  # Functions and types in generic.jl
  #----------------------------------
  export set_simulate_empirical_moments!, construct_objective_function!
  export set_priors!, set_empirical_moments!, set_bbSetup!, generate_bbSearchRange

  # Functions and types in save_load.jl
  #------------------------------------
  export read_priors, read_empirical_moments
  export saveSMMOptim, loadSMMOptim

  # Functions and types in optimize.jl
  #-----------------------------------
  export smmoptimize!, smm_minimizer




end # module
