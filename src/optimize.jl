"""
  smmoptimize!(sMMProblem::SMMProblem; verbose::Bool = true)

Function to launch an optimization. To be used after the following functions
have been called: (i) set_empirical_moments! (ii) set_priors!
(iii) set_simulate_empirical_moments! (iv) construct_objective_function!
"""

# [TODO] Instead of loading, just save the function periodically
function smmoptimize!(sMMProblem::SMMProblem; verbose::Bool = true)

  # Initialize a BlackBoxOptim problem
  # this modifies sMMProblem.bbSetup
  #-----------------------------------
  set_bbSetup!(sMMProblem)

  # Run the optimization by "batches"
  numberBatches = sMMProblem.options.maxFuncEvals/sMMProblem.options.saveSteps
  # Store best fitness and best candidates
  listBestFitness = []
  listBestCandidates = []
  @showprogress for i=1:numberBatches

      if verbose == true
        println("\n -------------------------------------")
        println("Batch $(i) / $(numberBatches)")
        println("------------------------------------- \n")
      end

      # at i=1, the optimization has not started yet
      # if i>1
      #   sMMProblem = loadSMMOptim(sMMProblem.options.saveName, verbose = verbose, verbose = true);
      # end

      # Run the optimization with BlackBoxOptim
      #----------------------------------------
      sMMProblem.bbResults = bboptimize(sMMProblem.bbSetup)

      # Save to disk
      #-------------
      saveSMMOptim(sMMProblem, verbose = verbose, saveName = sMMProblem.options.saveName, verbose = true);
      push!(listBestFitness, best_fitness(sMMProblem.bbResults))
      push!(listBestCandidates, best_candidate(sMMProblem.bbResults))
  end


  return listBestFitness, listBestCandidates


end
