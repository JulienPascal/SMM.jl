"""
  smmoptimize!(sMMProblem::SMMProblem; verbose::Bool = true)

Function to launch an optimization. To be used after the following functions
have been called: (i) set_empirical_moments! (ii) set_priors!
(iii) set_simulate_empirical_moments! (iv) construct_objective_function!
"""
function smmoptimize!(sMMProblem::SMMProblem; verbose::Bool = true)

  # Initialize a BlackBoxOptim problem
  # this modifies sMMProblem.bbSetup
  #-----------------------------------
  set_global_optimizer!(sMMProblem)

  # If the global optimizer is using BlackBoxOptim
  #-----------------------------------------------
  if is_bb_optimizer(sMMProblem.options.globalOptimizer) == true

    # Run the optimization by "batches"
    numberBatches = sMMProblem.options.maxFuncEvals/sMMProblem.options.saveSteps
    # Store best fitness and best candidates
    listBestFitness = []
    listBestCandidates = []
    @showprogress for i=1:numberBatches

        if verbose == true
          info("-------------------------------------")
          info("Batch $(i) / $(numberBatches)")
          info("------------------------------------- \n")
        end

        # Run the optimization with BlackBoxOptim
        #----------------------------------------
        sMMProblem.bbResults = bboptimize(sMMProblem.bbSetup)

        # Save to disk
        #-------------
        saveSMMOptim(sMMProblem, verbose = verbose, saveName = sMMProblem.options.saveName);
        push!(listBestFitness, best_fitness(sMMProblem.bbResults))
        push!(listBestCandidates, best_candidate(sMMProblem.bbResults))
    end

  # In the future, we may use other global minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end


  return listBestFitness, listBestCandidates


end

"""
  function smm_minimizer(sMMProblem::SMMProblem)

Function to get the parameter value minimizing the objective function
"""
function smm_minimizer(sMMProblem::SMMProblem)

  # If the global optimizer is using BlackBoxOptim
  #-----------------------------------------------
  if is_bb_optimizer(sMMProblem.options.globalOptimizer) == true

    best_candidate(sMMProblem.bbResults)

  # In the future, we may use other global minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end

end

"""
  smm_refine_globalmin!(sMMProblem::SMMProblem; verbose::Bool = true)

Function to refine the global minimum using a local minimization routine.
To be used after the following functions have been called: (i) set_empirical_moments!
(ii) set_priors! (iii) set_simulate_empirical_moments! (iv) construct_objective_function!
(v) smmoptimize!
"""
function smm_refine_globalmin!(sMMProblem::SMMProblem; verbose::Bool = true)


  x0 = smm_minimizer(sMMProblem)

  # Let's use the result from the global maximizer as the starting value
  #---------------------------------------------------------------------
  if verbose == true
    info("Refining the global maximum using a local algorithm.")
    info("Using Fminbox = $(sMMProblem.options.minBox)")
    info("Starting value = $(x0)")
  end


  if is_local_optimizer(sMMProblem.options.localOptimizer) == true

    # If using Fminbox option is true
    #-------------------------------
    if sMMProblem.options.minBox == true

      lower = create_lower_bound(sMMProblem)
      upper = create_upper_bound(sMMProblem)

      sMMProblem.optimResults = optimize(sMMProblem.objective_function, x0, lower, upper,
                                        iterations = sMMProblem.options.maxFuncEvals,
                                        convert_to_fminbox(sMMProblem.options.localOptimizer))

    else
      sMMProblem.optimResults = optimize(sMMProblem.objective_function, x0,
                                        convert_to_optim_algo(sMMProblem.options.localOptimizer),
                                        Optim.Options(iterations = sMMProblem.options.maxFuncEvals))
    end

  # In the future, we may use other local minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.localOptimizer = $(sMMProblem.options.localOptimizer) is not supported.")

  end

  return Optim.minimizer(sMMProblem.optimResults)

end


"""
  function smm_minimizer(sMMProblem::SMMProblem)

Function to get the parameter value minimizing the objective function (local)
"""
function smm_local_minimizer(sMMProblem::SMMProblem)

  # If the global optimizer is using BlackBoxOptim
  #-----------------------------------------------
  if is_optim_optimizer(sMMProblem.options.localOptimizer) == true

    Optim.minimizer(sMMProblem.optimResults)

  # In the future, we may use other global minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end

end

"""
  function smm_minimum(sMMProblem::SMMProblem)

Function to get the local minimum value of the objetive function
"""
function smm_local_minimum(sMMProblem::SMMProblem)

  # If the global optimizer is using BlackBoxOptim
  #-----------------------------------------------
  if is_optim_optimizer(sMMProblem.options.localOptimizer) == true

    Optim.minimum(sMMProblem.optimResults)

  # In the future, we may use other global minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end

end


"""
  local_to_global!(sMMProblem::SMMProblem; x0 = Array{Float64}(0,0), nums::Int64 = nworkers(), verbose::Bool = true)

Function to run several local minimization algorithms in parallel, with different
starting values. The minimum is calculated as the minimum of the local minima.
Changes sMMProblem.optimResults. This function also returns
a list containing Optim results.
"""
function local_to_global!(sMMProblem::SMMProblem; x0 = Array{Float64}(0,0), nums::Int64 = nworkers(), verbose::Bool = true)

  # Safety checks
  #--------------
  if nums < nworkers()
    Base.errors("nums < nworkers()")
  elseif nums > nworkers()
    info("nums > nworkers(). Some starting values will be ignored.")
  end

  # To store minization results
  #----------------------------
  results = []

  # Look for valid starting values (for which convergence is reached)
  #-------------------------------------------------------------------
  if x0 == Array{Float64}(0,0)
    myGrid = search_starting_values(sMMProblem, nums, verbose = verbose)
  # Using starting values provided by the user
  #-------------------------------------------
  else

    # Check that enough starting values were provided
    if size(x0, 1) < nworkers()
      info("$(size(x0, 1)) starting value(s) were provided.")
      Base.error("The minimum number of starting value(s) to provide is $(nworkers()).")
    end

    myGrid = x0

  end


  # If the local optimizer is using Optim
  #--------------------------------------
  if is_optim_optimizer(sMMProblem.options.localOptimizer) == true

      # A. Starting tasks on available workers
      #---------------------------------------
      @sync for (workerIndex, w) in enumerate(workers())

        @async push!(results, @fetchfrom w wrap_smm_localmin(sMMProblem, myGrid[workerIndex,:], verbose = true))

      end


    # B. Looking for the minimum
    #----------------------------
    # Initialization
    minIndex = 0
    minValue = Inf
    minimizerValue = zeros(length(keys(sMMProblem.priors)))
    nbConvergenceReached = 0
    listOptimResults = []

    for (workerIndex, w) in enumerate(workers())

      push!(listOptimResults, results[workerIndex])

      try

        minimumValue = Optim.minimum(results[workerIndex])
        minimizer = Optim.minimizer(results[workerIndex])

        if minimumValue < minValue && Optim.converged(results[workerIndex]) == true

          minIndex = workerIndex
          minValue = minimumValue
          minimizerValue = minimizer
          nbConvergenceReached += 1

        end

      catch myError
        info("$(myError)")
      end

    end


    # D. If none of the optimization converged
    #-----------------------------------------
    if minIndex == 0
      info("None of the local optimizer algorithm converged.")
    else
      info("Convergence reached for $(nbConvergenceReached) worker(s).")
      info("Minimum value found with worker $(minIndex)")
      sMMProblem.optimResults = results[minIndex]
    end

  # In the future, we may use other local minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.localOptimizer = $(sMMProblem.options.localOptimizer) is not supported.")

  end

  if verbose == true
    if nbConvergenceReached != 0
      info("Best value found with starting values = $(myGrid[minIndex,:]).")
      info("Best value = $(minValue).")
      info("Minimizer = $(minimizerValue)")
    end
  end

  return listOptimResults

end


"""
  smm_localmin(sMMProblem::SMMProblem, x0::Array{Float64,1}; verbose::Bool = true)

Function find a local minimum using a local minimization routine, with starting value x0.
To be used after the following functions have been called: (i) set_empirical_moments!
(ii) set_priors! (iii) set_simulate_empirical_moments! (iv) construct_objective_function!
"""
function smm_localmin(sMMProblem::SMMProblem, x0::Array{Float64,1}; verbose::Bool = true)

  # Let's use the result from the global maximizer as the starting value
  #---------------------------------------------------------------------
  if verbose == true
    info("Starting value = $(x0)")
    info("Using Fminbox = $(sMMProblem.options.minBox)")
  end


  if is_local_optimizer(sMMProblem.options.localOptimizer) == true

    # If using Fminbox option is true
    #-------------------------------
    if sMMProblem.options.minBox == true

      lower = create_lower_bound(sMMProblem)
      upper = create_upper_bound(sMMProblem)

      optimResults = optimize(sMMProblem.objective_function, x0, lower, upper,
                                        iterations = sMMProblem.options.maxFuncEvals,
                                        convert_to_fminbox(sMMProblem.options.localOptimizer))
    else
      optimResults = optimize(sMMProblem.objective_function, x0,
                            convert_to_optim_algo(sMMProblem.options.localOptimizer),
                            Optim.Options(iterations = sMMProblem.options.maxFuncEvals))
    end

  # In the future, we may use other local minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.localOptimizer = $(sMMProblem.options.localOptimizer) is not supported.")

  end

  return optimResults

end


"""

"""
function wrap_smm_localmin(sMMProblem::SMMProblem, x0::Array{Float64,1}; verbose::Bool = true)

  try
    smm_localmin(sMMProblem, x0, verbose = verbose)
  catch myError
    info("$(myError)")
    info("Error with smm_localmin")
  end

end


"""
  search_starting_values(sMMProblem::SMMProblem; verbose::Bool = true, a::Array{Float64,1}, b::Array{Float64,1}, nums::Int64)

Search for nums valid starting values. To be used after the following functions have been called:
(i) set_empirical_moments! (ii) set_priors! (iii) set_simulate_empirical_moments!
(iv) construct_objective_function!
"""
function search_starting_values(sMMProblem::SMMProblem, numPoints::Int64; verbose::Bool = true, maxNbTrials::Int64 = 1000)

  # Safety Check
  #-------------
  if is_optim_optimizer(sMMProblem.options.localOptimizer) == false
    Base.error("sMMProblem.options.localOptimizer = $(sMMProblem.options.localOptimizer) is not supported.")
  end

  if verbose == true
    info("Searching for $(numPoints) valid starting values")
  end

  # Generate upper and lower bounds vector readable by create_grid_stochastic
  #--------------------------------------------------------------------------
  lower_bound = zeros(length(keys(sMMProblem.priors)))
  upper_bound = zeros(length(keys(sMMProblem.priors)))

  for (kIndex, k) in enumerate(keys(sMMProblem.priors))
    lower_bound[kIndex] = sMMProblem.priors[k][2]
    upper_bound[kIndex] = sMMProblem.priors[k][3]
  end

  #Each row is a new point and each column is a dimension of this points.
  #---------------------------------------------------------------------
  Validx0 = zeros(numPoints, length(lower_bound))
  nbValidx0Found = 0

  # Create many grids (stochastic draws) with many potential points
  #----------------------------------------------------------------
  if verbose == true
    info("Creating $(maxNbTrials) potential grids")
  end

  listGrids = []
  for i=1:maxNbTrials
    push!(listGrids, create_grid_stochastic(lower_bound, upper_bound, numPoints, gridType = :uniform))
  end

  listGridsIndex = 0

  # Looping until all the points have been found
  #---------------------------------------------
  while nbValidx0Found < numPoints

    results = []
    listGridsIndex += 1

    if listGridsIndex > maxNbTrials
      Base.error("Maximum number of attempts reached without success. maxNbTrials = $(maxNbTrials)")
    end

    # Use available workers to simulate moments
    #------------------------------------------
    @sync for (workerIndex, w) in enumerate(workers())

      @async push!(results, @fetchfrom w sMMProblem.objective_function(listGrids[listGridsIndex][workerIndex, :]))

    end

    # Check for convergence
    #----------------------
    for (workerIndex, w) in enumerate(workers())

      try

        distanceValue = results[workerIndex]

        if isinf(distanceValue) == false && distanceValue != 999999.0

          nbValidx0Found +=1

          if nbValidx0Found <= numPoints

            Validx0[nbValidx0Found,:] = listGrids[listGridsIndex][workerIndex, :]
            info("Valid starting value = $(Validx0[nbValidx0Found,:])")
          end

        end

      catch myError
        info("$(myError)")
      end

    end

  end

  if verbose == true
    info("Found $(nbValidx0Found) valid starting values")
  end

  return Validx0

end
