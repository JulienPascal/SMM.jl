
"""
# SMM.jl
This is the documentation for SMM.jl

"""

module SMM

  # Dependencies
  #-------------
  using BlackBoxOptim
  using Optim
  using JLD2
  using Plots
  using CSV
  using ProgressMeter
  using Calculus
  using CompEcon
  using Distributions
  using DataFrames

  @static if VERSION < v"0.7.0"
      using DataStructures
  else
      using OrderedCollections
  end



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

  # Functions to do inference
  #---------------------------------------------
  include("econometrics.jl")

  # Functions to do plots
  #---------------------------------------------
  include("plots.jl")

  # Exports
  #--------
  # Functions and types in types.jl
  #----------------------------------
  export SMMOptions, SMMProblem
  export default_function, rosenbrock2d
  export is_global_optimizer, is_local_optimizer
  export convert_to_optim_algo, convert_to_fminbox
  export is_bb_optimizer, is_optim_optimizer

  # Functions and types in api.jl
  #-------------------------------

  # Functions and types in generic.jl
  #----------------------------------
  export set_simulate_empirical_moments!, construct_objective_function!
  export set_priors!, set_empirical_moments!, set_Sigma0!
  export set_bbSetup!, generate_bbSearchRange
  export create_lower_bound, create_upper_bound
  export set_global_optimizer!
  export create_grid, create_grid_stochastic, generate_std, latin_hypercube_sampling

  # Functions and types in save_load.jl
  #------------------------------------
  export read_priors, read_empirical_moments
  export saveSMMOptim, loadSMMOptim

  # Functions and types in optimize.jl
  #-----------------------------------
  export smmoptimize!, smm_minimizer
  export smm_refine_globalmin!, smm_local_minimizer
  export smm_local_minimum
  export smm_localmin, local_to_global!

  # Functions in econometrics.jl
  #-----------------------------
  export calculate_D, calculate_Avar!, calculate_se, calculate_t, calculate_pvalue, calculate_CI
  export summary_table

  # Functions in plots.jl
  #----------------------
  export smm_slices, drawMe




end # module
