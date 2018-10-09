struct SMMOptions
	bbOptimizer::Symbol
end

function SMMOptions(;bbOptimizer::Symbol=:adaptive_de_rand_1_bin_radiuslimited)

	SMMOptions(bbOptimizer)
end

mutable struct SMMProblem
	iter::Int64
end
