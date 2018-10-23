"""
	SMMOptions

SMMOptions is a struct that contains options related to the optimization

"""
struct SMMOptions
	bbOptimizer::Symbol
	maxFuncEvals::Int64
	saveSteps::Int64
	saveName::String
end

function SMMOptions( ;bbOptimizer::Symbol=:adaptive_de_rand_1_bin_radiuslimited,
											maxFuncEvals::Int64=1000,
											saveSteps::Int64 = 100,
											saveName::String = string(Dates.now()))

	# Safety Checks
	#--------------
	if saveSteps > maxFuncEvals
		error("Error in the constructor for SMMOptions. \n saveSteps = $(saveSteps) > maxFuncEvals = $(maxFuncEvals)")
	end

	if mod(maxFuncEvals, saveSteps) != 0
		error("Error in the constructor for SMMOptions. \n maxFuncEvals should be a multiple of saveSteps")
	end


	SMMOptions(bbOptimizer,
						maxFuncEvals,
						saveSteps,
						saveName)

end

"""
	SMMProblem

SMMProblem is a mutable struct that caries all the information needed to
perform the optimization and display the results

"""
mutable struct SMMProblem
	iter::Int64
	priors::OrderedDict{String,Array{Float64,1}}
	empiricalMoments::OrderedDict{String,Array{Float64,1}}
	simulatedMoments::OrderedDict{String, Float64}
	distanceEmpSimMoments::Float64
	simulate_empirical_moments::Function
	objective_function::Function
	options::SMMOptions
	bbSetup::BlackBoxOptim.OptController
	bbResults::BlackBoxOptim.OptimizationResults
end

# Constructor for SMMProblem
#------------------------------------------------------------------------------
function SMMProblem(  ;iter::Int64 = 0,
											priors::OrderedDict{String,Array{Float64,1}} = OrderedDict{String,Array{Float64,1}}(),
											empiricalMoments::OrderedDict{String,Array{Float64,1}} = OrderedDict{String,Array{Float64,1}}(),
											simulatedMoments::OrderedDict{String, Float64} = OrderedDict{String,Float64}(),
											distanceEmpSimMoments::Float64 = 0.,
											simulate_empirical_moments::Function = default_function,
											objective_function::Function = default_function,
											options::SMMOptions = SMMOptions(),
											bbSetup::BlackBoxOptim.OptController = defaultbbOptimOptController,
											bbResults::BlackBoxOptim.OptimizationResults = defaultbbOptimOptimizationResults)

	SMMProblem(iter,
						priors,
						empiricalMoments,
						simulatedMoments,
						distanceEmpSimMoments,
						simulate_empirical_moments,
						objective_function,
						options,
						bbSetup,
						bbResults)

end

"""
	default_function(x)

Function x->x. Used to initialize functions
"""
function default_function(x)
	println("default_function, returns input")
	x
end


"""
	rosenbrock2d(x)

Rosenbrock function. Used to initialize BlackBoxOptim.OptController().
"""
function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

# It is quite useful to have "default" BlackBoxOptim.OptController and
# BlackBoxOptim.OptimizationResults objects, since
# BlackBoxOptim.OptController() and BlackBoxOptim.OptimizationResults()
# do not work
#-------------------------------------------------------------------------------
defaultbbOptimOptController = bbsetup(x -> rosenbrock2d(x);
																		Method=:adaptive_de_rand_1_bin_radiuslimited,
																		SearchRange = (-5.0, 5.0),
             												NumDimensions = 2, MaxFuncEvals = 2,
																		TraceMode = :silent)

defaultbbOptimOptimizationResults = bboptimize(defaultbbOptimOptController)
