
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


  # Exports
  #--------
  # Functions and types in "types.jl"
  #----------------------------------
  export SMMOptions, SMMProblem
  export default_function, rosenbrock2d

  # Functions and types in "api.jl"
  #-------------------------------

  # Functions and types in generic.jl
  #-----------------------------------------
  export set_simulate_empirical_moments!, construct_objective_function!




end # module
