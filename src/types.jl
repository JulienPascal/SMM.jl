"""
	SMMOptions

SMMOptions is a struct that contains options related to the optimization.
"""
struct SMMOptions
	globalOptimizer::Symbol #algorithm for finding a global maximum
	localOptimizer::Symbol 	#algorithm for finding a local maximum
	maxFuncEvals::Int64			#maximum number of evaluations
	saveSteps::Int64				#maximum number of steps
	saveName::String				#name under which the optimization should be saved
	showDistance::Bool			#show the distance, everytime the objective function is calculated?
end

function SMMOptions( ;globalOptimizer::Symbol=:dxnes,
											localOptimizer::Symbol=:LBFGS,
											maxFuncEvals::Int64=1000,
											saveSteps::Int64 = 100,
											saveName::String = string(Dates.now()),
											showDistance::Bool = false)

	# Safety Checks
	#--------------
	if saveSteps == 0
		error("You cannot set saveSteps equal to 0.")
	end

	if saveSteps > maxFuncEvals
		error("Error in the constructor for SMMOptions. \n saveSteps = $(saveSteps) > maxFuncEvals = $(maxFuncEvals)")
	end

	if mod(maxFuncEvals, saveSteps) != 0
		error("Error in the constructor for SMMOptions. \n maxFuncEvals should be a multiple of saveSteps")
	end

	# Warnings
	#----------
	# Saving too often leads to a poor performance if using BlackBoxOptim
	# see here: https://github.com/robertfeldt/BlackBoxOptim.jl/issues/99
	#--------------------------------------------------------------------
	if maxFuncEvals/saveSteps < 6 && is_bb_optimizer(globalOptimizer) == true
		info("WARNING. When maxFuncEvals/saveSteps < 6 using BlackBoxOptim, the performance
						of the global maximizer may deteriorate.")
	end




	SMMOptions(globalOptimizer,
						localOptimizer,
						maxFuncEvals,
						saveSteps,
						saveName,
						showDistance)

end

"""
	SMMProblem

SMMProblem is a mutable struct that caries all the information needed to
perform the optimization and display the results.
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
	bbSetup::BlackBoxOptim.OptController					#set up when using BlackBoxOptim (global minimum)
	bbResults::BlackBoxOptim.OptimizationResults	#results when using BlackBoxOptim (global minimum)
	optimResults::Optim.OptimizationResults				#results when using Optim (local minimum)
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
											bbResults::BlackBoxOptim.OptimizationResults = defaultbbOptimOptimizationResults,
											optimResults::Optim.OptimizationResults = defaultOptimResults)

	SMMProblem(iter,
						priors,
						empiricalMoments,
						simulatedMoments,
						distanceEmpSimMoments,
						simulate_empirical_moments,
						objective_function,
						options,
						bbSetup,
						bbResults,
						optimResults)

end

"""
	default_function(x)

Function x->x. Used to initialize functions.
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
																		Method=:dxnes,
																		SearchRange = (-5.0, 5.0),
             												NumDimensions = 2, MaxFuncEvals = 2,
																		TraceMode = :silent)

defaultbbOptimOptimizationResults = bboptimize(defaultbbOptimOptController)

defaultOptimResults = optimize(rosenbrock2d, [0.0, 0.0], LBFGS())

"""
	is_global_optimizer(s::Symbol)

function to check that the global optimizer chosen is available.
"""
function is_global_optimizer(s::Symbol)

	# Global Optimizers using BlackBoxOptim
	# source: https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/examples/benchmarking/latest_toplist.csv
	#------------------------------------------------------------------------------
	listValidGlobalOptimizers = [:dxnes, :adaptive_de_rand_1_bin_radiuslimited, :xnes,
															 :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin,
															 :generating_set_search, :de_rand_1_bin,
															 :separable_nes, :resampling_inheritance_memetic_search,
															 :probabilistic_descent, :resampling_memetic_search,
															 :de_rand_2_bin_radiuslimited, :de_rand_2_bin,
															 :random_search, :simultaneous_perturbation_stochastic_approximation]

	in(s, listValidGlobalOptimizers)

end

"""
	is_bb_optimizer(s::Symbol)

function to check whether the global optimizer is using BlackBoxOptim
"""
function is_bb_optimizer(s::Symbol)

	# source: https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/examples/benchmarking/latest_toplist.csv
	listbbOptimizers = [:dxnes, :adaptive_de_rand_1_bin_radiuslimited, :xnes,
															 :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin,
															 :generating_set_search, :de_rand_1_bin,
															 :separable_nes, :resampling_inheritance_memetic_search,
															 :probabilistic_descent, :resampling_memetic_search,
															 :de_rand_2_bin_radiuslimited, :de_rand_2_bin,
															 :random_search, :simultaneous_perturbation_stochastic_approximation]

	in(s, listbbOptimizers)

end

"""
	is_local_optimizer(s::Symbol)

function to check that the local optimizer chosen is available.
"""
function is_local_optimizer(s::Symbol)

	# source: https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/examples/benchmarking/latest_toplist.csv
	listValidLocalOptimizers = [:NelderMead, :SimulatedAnnealing, :ParticleSwarm,
															:BFGS, :LBFGS, :ConjugateGradient, :GradientDescent,
															:MomentumGradientDescent, :AcceleratedGradientDescent]

	in(s, listValidLocalOptimizers)

end

"""
	is_optim_optimizer(s::Symbol)

function to check whether the local minimizer uses the package Optim.
"""
function is_optim_optimizer(s::Symbol)

	# source: https://github.com/robertfeldt/BlackBoxOptim.jl/blob/master/examples/benchmarking/latest_toplist.csv
	listOptimOptimizers = [:NelderMead, :SimulatedAnnealing, :ParticleSwarm,
															:BFGS, :LBFGS, :ConjugateGradient, :GradientDescent,
															:MomentumGradientDescent, :AcceleratedGradientDescent]

	in(s, listOptimOptimizers)

end


"""
	convert_to_optim_algo(s::Symbol)

function to convert local optimizer (of type Symbol) to an Optim algo.
"""
function convert_to_optim_algo(s::Symbol)

	if s == :NelderMead

		output = NelderMead()

	elseif s == :SimulatedAnnealing

		output = SimulatedAnnealing()

	elseif s == :ParticleSwarm

		output = ParticleSwarm()

	elseif s == :BFGS

		output = BFGS()

	elseif s == :LBFGS

		output = LBFGS()

	elseif s == :ConjugateGradient

		output = ConjugateGradient()

	elseif s == :GradientDescent

		output = GradientDescent()

	elseif s == :MomentumGradientDescent

		output = MomentumGradientDescent()

	elseif s == :AcceleratedGradientDescent

		output = AcceleratedGradientDescent()

	else

		Base.error("$(s) is not in the list of algorithm supported by Optim.")

	end

	return output

end
