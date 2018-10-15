"""
  set_simulate_empirical_moments!(sMMProblem::SMMProblem, f::Function)

Function to set the field simulate_empirical_moments for a SMMProblem.
The function simulate_empirical_moments takes parameter values and return
the corresponding simulate moments values.

"""
function set_simulate_empirical_moments!(sMMProblem::SMMProblem, f::Function)

  sMMProblem.simulate_empirical_moments = f

end

"""
  construct_objective_function!(sMMProblem::SMMProblem, f::Function)

Function that construct an objective function, using the function
SMMProblem.simulate_empirical_moments. For the moment, the objective function
returns the mean of the percentage deviation of simulated moments from
their true value.

"""
function construct_objective_function!(sMMProblem::SMMProblem, objectiveType::Symbol = :percent)


    function objective_function_percent(x)


      # Initialization
      #----------------
      distanceEmpSimMoments = 999999.0

      #------------------------------------------------------------------------
      # A.
      # Generate simulated moments
      #------------------------------------------------------------------------
      simulatedMoments, convergence = try

          sMMProblem.simulate_empirical_moments(x), 1

      catch errorSimulation

        info("An error occured with parameter values = $(x)")
        info("$(errorSimulation)")

          OrderedDict{String,Array{Float64,1}}(), 0

      end

      #------------------------------------------------------------------------
      # B.
      # If generating moment was successful, calculate distance between empirical
      # and simulated moments
      # (If no convergence, returns large penalty value 999999.0)
      #------------------------------------------------------------------------
      if convergence == 1

        # to store the distance between empirical and simulated moments
        arrayDistance = zeros(length(keys(sMMProblem.empiricalMoments)))

        for (indexMoment, k) in enumerate(keys(sMMProblem.empiricalMoments))

          # * sMMProblem.empiricalMoments[k][1] is the empirical moments
          # * sMMProblem.empiricalMoments[k][2] is the weight associated to this moment
          #---------------------------------------------------------------------
          arrayDistance[indexMoment] = ((sMMProblem.empiricalMoments[k][1] - simulatedMoments[k]) /sMMProblem.empiricalMoments[k][2])^2

        end

        distanceEmpSimMoments = mean(arrayDistance)

      end

      return distanceEmpSimMoments

    end


  if objectiveType == :percent
    sMMProblem.objective_function = objective_function_percent
  else
    error("Error in construct_objective_function. objectiveType currently only support :percent.")
  end


end

"""
  set_priors!(sMMProblem::SMMProblem, priors::DataStructures.OrderedDict{String,Array{Float64,1}})

Function to change the field sMMProblem.priors

"""
function set_priors!(sMMProblem::SMMProblem, priors::DataStructures.OrderedDict{String,Array{Float64,1}})

  sMMProblem.priors = priors

end

"""
   set_empirical_moments!(sMMProblem::SMMProblem, empiricalMoments::DataStructures.OrderedDict{String,Array{Float64,1}})

Function to change the field sMMProblem.empiricalMoments

"""
function set_empirical_moments!(sMMProblem::SMMProblem, empiricalMoments::DataStructures.OrderedDict{String,Array{Float64,1}})

  sMMProblem.empiricalMoments = empiricalMoments

end

"""
  set_bbSetup!(sMMProblem::SMMProblem)

Function to set the field bbSetup for a SMMProblem.
"""
function set_bbSetup!(sMMProblem::SMMProblem)

  # A. using sMMProblem.priors, generate searchRange:
  #-------------------------------------------------
  mySearchRange = generate_bbSearchRange(sMMProblem)

  info("$(nworkers()) worker(s) detected")

  if nworkers() == 1
    info("Starting optimization in serial")
    sMMProblem.bbSetup = bbsetup(sMMProblem.objective_function;
                              Method = sMMProblem.options.bbOptimizer,
                              SearchRange = mySearchRange,
                              MaxFuncEvals = sMMProblem.options.saveSteps,
                              TraceMode = :verbose,
                              NumDimensions = length(keys(sMMProblem.priors)))
  else
    info("Starting optimization in parallel")
    sMMProblem.bbSetup = bbsetup(sMMProblem.objective_function;
                                Method = sMMProblem.options.bbOptimizer,
                                SearchRange = mySearchRange,
                                MaxFuncEvals = sMMProblem.options.saveSteps,
                                Workers = workers(),
                                TraceMode = :verbose,
                                NumDimensions = length(keys(sMMProblem.priors)))
  end


end

"""
  generate_bbSearchRange(sMMProblem::SMMProblem)

Function to generate a search range that matches the convention used by
BlackBoxOptim.
"""
function generate_bbSearchRange(sMMProblem::SMMProblem)

  # sMMProblem.priors["key"][1] contains the initial guess
  # sMMProblem.priors["key"][2] contains the lower bound
  # sMMProblem.priors["key"][3] contains the upper bound
  #-----------------------------------------------------
  [(sMMProblem.priors[k][2], sMMProblem.priors[k][3]) for k in keys(sMMProblem.priors)]
end
