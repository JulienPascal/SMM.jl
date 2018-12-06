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
    info("Starting value = $(x0)")
  end


  if is_local_optimizer(sMMProblem.options.localOptimizer) == true

    sMMProblem.optimResults = optimize(sMMProblem.objective_function, x0,
                                      convert_to_optim_algo(sMMProblem.options.localOptimizer),
                                      Optim.Options(iterations = sMMProblem.options.maxFuncEvals))

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


#=
"""
  local_to_global(sMMProblem::SMMProblem)

Function to run
"""
function local_to_global(sMMProblem::SMMProblem; nums::Int64 = 2, verbose::Bool = true)

  # To store minization results
  #----------------------------
  results = []
  refs = []
  refIndex = 1


  # Create a grid for starting points
  #----------------------------------



  # If the global optimizer is using BlackBoxOptim
  #-----------------------------------------------
  if is_optim_optimizer(sMMProblem.options.localOptimizer) == true

    # A. Spawning
    #------------
    @sync for w in workers()

      @async refs[refIndex] = @spawnat w

    end

    # B. Fetching
    #-------------
    for w in workers()

      push!(results, @fetch refs[refIndex])

    end

    # C. Looking for the minimum
    #----------------------------
    # Initialization
    minIndex = 0
    minValue = Inf

    for w in workers()

      minimumValue = Optim.minimum(results[w])

      if minimumValue < minValue && Optim.converged(results[w]) == true

        minIndex = w
        minValue = minimumValue

      end

    end

    # D. If none of the optimization converged
    #-----------------------------------------
    if minIndex == 0
      info("None of the local optimizer algorithm conerged.")
    end

  # In the future, we may use other global minimizer
  # routines. For the moment, let's return an error
  #-------------------------------------------------
  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end

end
=#
