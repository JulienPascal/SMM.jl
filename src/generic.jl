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
  set_priors!(sMMProblem::SMMProblem, priors::OrderedDict{String,Array{Float64,1}})

Function to change the field sMMProblem.priors

"""
function set_priors!(sMMProblem::SMMProblem, priors::OrderedDict{String,Array{Float64,1}})

  sMMProblem.priors = priors

end

"""
   set_empirical_moments!(sMMProblem::SMMProblem, empiricalMoments::OrderedDict{String,Array{Float64,1}})

Function to change the field sMMProblem.empiricalMoments

"""
function set_empirical_moments!(sMMProblem::SMMProblem, empiricalMoments::OrderedDict{String,Array{Float64,1}})

  sMMProblem.empiricalMoments = empiricalMoments

end

"""
  set_global_optimizer!(sMMProblem::SMMProblem)

Function to set the fields corresponding to the global
optimizer problem.
"""
function set_global_optimizer!(sMMProblem::SMMProblem)

  if is_bb_optimizer(sMMProblem.options.globalOptimizer) == true

    set_bbSetup!(sMMProblem)

  else

    Base.error("sMMProblem.options.globalOptimizer = $(sMMProblem.options.globalOptimizer) is not supported.")

  end

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
   # Debug:
   #-------
   @sync for (idx, pid) in enumerate(workers())
     @async @spawnat(pid, println("hello"))
   end

  if nworkers() == 1
    info("Starting optimization in serial")
    sMMProblem.bbSetup = bbsetup(sMMProblem.objective_function;
                              Method = sMMProblem.options.globalOptimizer,
                              SearchRange = mySearchRange,
                              MaxFuncEvals = sMMProblem.options.saveSteps,
                              TraceMode = :verbose,
                              NumDimensions = length(keys(sMMProblem.priors)))
  else
    info("Starting optimization in parallel")
    sMMProblem.bbSetup = bbsetup(sMMProblem.objective_function;
                                Method = sMMProblem.options.globalOptimizer,
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


"""
  cartesian_grid(a::Array{Float64,1}, b::Array{Float64,1}, nums::Int64)

Function to create a regular cartesian grid. a is a vector of lower
bounds, b a vector of upper bounds and nums is the number of points
along each dimension. This function works for up to 15 dimensions.
Returns a NTuple.
"""
function cartesian_grid(a::Array{Float64,1}, b::Array{Float64,1}, nums::Int64)

  #Safety checks
  #-------------
  if nums < 2
    Base.error("The input nums should be >= 2. nums = $(nums).")
  end

  if length(a) != length(b)
    Base.error("length(a) != length(b)")
  end

  nodes = [collect(linspace(a[i], b[i], nums)) for i in 1:length(a)]

  # There is probably a better way of doing this:
  #----------------------------------------------
  if length(a) == 1
    points = collect(Iterators.product(nodes[1]))
  elseif length(a) == 2
    points = collect(Iterators.product(nodes[1], nodes[2]))
  elseif length(a) == 3
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3]))
  elseif length(a) == 4
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4]))
  elseif length(a) == 5
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5]))
  elseif length(a) == 6
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6]))
  elseif length(a) == 7
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7]))
  elseif length(a) == 8
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8]))
  elseif length(a) == 9
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9]))
  elseif length(a) == 10
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10]))
  elseif length(a) == 11
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10],
                      nodes[11]))
  elseif length(a) == 12
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10],
                      nodes[11], nodes[12]))
  elseif length(a) == 13
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10],
                      nodes[11], nodes[12], nodes[13]))
  elseif length(a) == 14
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10],
                      nodes[11], nodes[12], nodes[13], nodes[14]))
  elseif length(a) == 15
    points = collect(Iterators.product(nodes[1], nodes[2], nodes[3], nodes[4],
                      nodes[5], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10],
                      nodes[11], nodes[12], nodes[13], nodes[14], nodes[15]))
  else

    Base.error("cartesian_grid can take up to 15 dimensions")

  end

  return points

end

"""
  create_grid(a::Array{Float64,1}, b::Array{Float64,1}, nums::Int64)

Function to create a grid. a is a vector of lower
bounds, b a vector of upper bounds and nums is the number of points
along each dimension. The type of grid to use can be specified using gridType.
By default, the functions returns the cartesian product. The output is an
Array{Float64,2}, where each row is a new point and each column is a dimension
of this points.
"""
function create_grid(a::Array{Float64,1}, b::Array{Float64,1}, nums::Int64; gridType::Symbol = :lin)

	  #Safety checks
	  #-------------
	  if nums < 2
	    Base.error("The input nums should be >= 2. nums = $(nums).")
	  end

	  if length(a) != length(b)
	    Base.error("length(a) != length(b)")
	  end

	if gridType != :cheb && gridType == :spli && gridType == :lin
		Base.error("gridType has to be either :chebn, :spli or :lin")
	end

  Fspace = fundefn(gridType, collect([nums for i =1:length(a)]), a, b)

	Fnodes = funnode(Fspace)[1]

end
