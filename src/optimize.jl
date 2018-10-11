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
  set_bbSetup!(sMMProblem)

  # Run the optimization by "batches"
  numberBatches = sMMProblem.options.maxFuncEvals/sMMProblem.options.saveSteps

  @showprogress for i=1:numberBatches

      if verbose == true
        println("\n -------------------------------------")
        println("Batch $(i) / $(numberBatches)")
        println("------------------------------------- \n")
      end

      # at i=1, the optimization has not started yet
      if i>1
        sMMProblem = loadSMMOptim(sMMProblem.options.saveName, verbose = verbose);
      end

      # Run the optimization with BlackBoxOptim
      #----------------------------------------
      sMMProblem.bbResults = bboptimize(sMMProblem.bbSetup)

      # Save to disk
      #-------------
      saveSMMOptim(sMMProblem, verbose = verbose, saveName = sMMProblem.options.saveName);

  end


  return sMMProblem


end
