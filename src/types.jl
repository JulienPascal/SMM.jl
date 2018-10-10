"""
	SMMOptions

SMMOptions is a struct that contains options related to the optimization

"""
struct SMMOptions
	bbOptimizer::Symbol
end

function SMMOptions(;bbOptimizer::Symbol=:adaptive_de_rand_1_bin_radiuslimited)

	SMMOptions(bbOptimizer)

end

"""
	SMMProblem

SMMProblem is a mutable struct that caries all the information needed to
perform the optimization and display the results

"""
mutable struct SMMProblem
	iter::Int64
	priors::Dict{String, Any}
	empiricalMoments::Dict{String, Any}
	simulatedMoments::Dict{String, Any}
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
											priors::Dict{String, Any} = Dict{String, Any}(),
											empiricalMoments::Dict{String, Any} = Dict{String, Any}(),
											simulatedMoments::Dict{String, Any} = Dict{String, Any}(),
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
	info("default_function, returns input")
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
             												NumDimensions = 2, MaxFuncEvals = 2);

defaultbbOptimOptimizationResults = bboptimize(defaultbbOptimOptController);
