@everywhere using SMM
@everywhere using Distributions
@everywhere using DataStructures

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# Tests that are too long to be sent to TRAVIS
#---------------------------------------------
# 6d problem
#--------------------------------------------------------------
@testset "Testing smmoptimize with 6d and same magnitude" begin

    # Rermark:
    # It is important NOT to use srand()
    # within the function simulate_empirical_moments!
    # Otherwise, BlackBoxOptim does not find the solution
    #----------------------------------------------------
    # The difference of magniture make it more difficult to find the minimum
    tol6dMean = 0.5

    @everywhere function functionTest6d(x)

        d = MvNormal(x, eye(6))

        output = OrderedDict{String,Float64}()
        
        draws = rand(d, 1000000)
        for i in 1:6
            output["mean$(i)"] = mean(draws[i,:])
        end

        return output
    end


    t = SMMProblem(options = SMMOptions(maxFuncEvals=1000,saveSteps = 1000, saveName = "parallel6D"))


    #------------------------------------------------------
    dictEmpiricalMoments = OrderedDict{String,Array{Float64,1}}()
    dictEmpiricalMoments["mean1"] = [ -1.0; -1.0]
    dictEmpiricalMoments["mean2"] = [ 1.0; 1.0]
    dictEmpiricalMoments["mean3"] = [ 0.5; 0.5]
    dictEmpiricalMoments["mean4"] = [ -0.5; -0.5]
    dictEmpiricalMoments["mean5"] = [ 0.7; 0.7]
    dictEmpiricalMoments["mean6"] = [-0.7; -0.7]
    set_empirical_moments!(t, dictEmpiricalMoments)

    dictPriors = OrderedDict{String,Array{Float64,1}}()
    dictPriors["mu1"] = [0.2,-3., 3.]
    dictPriors["mu2"] = [-0.2,-2.,2.]
    dictPriors["mu3"] = [-0.3,-2.,2.]
    dictPriors["mu4"] = [-0.4,-2., 2.]
    dictPriors["mu5"] = [0.3,-2., 2.]
    dictPriors["mu6"] = [0.4,-2., 2.]

    set_priors!(t, dictPriors)

    # A. Set the function: parameter -> simulated moments
    set_simulate_empirical_moments!(t, functionTest6d)

    # B. Construct the objective function, using the function: parameter -> simulated moments
    # and moments' weights:
    construct_objective_function!(t)

    # C. Run the optimization
    # This function first modifies t.bbSetup
    # and then modifies t.bbResults

    # [TODO] loading and restarting does
    # not seem to work when using workers in parallel
    # Figure out why
    listBestFitness, listBestCandidates = smmoptimize!(t, verbose = true)

    @test best_candidate(t.bbResults)[1] ≈ -1.0 atol = tol6dMean
    @test best_candidate(t.bbResults)[2] ≈ 1.0 atol = tol6dMean
    @test best_candidate(t.bbResults)[3] ≈ 0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[4] ≈ -0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[5] ≈ 0.7 atol = tol6dMean
    @test best_candidate(t.bbResults)[6] ≈ -0.7 atol = tol6dMean

end


# Ambitious problem, with 18 parameters
# different magnitudes
#--------------------------------------
# source: https://github.com/floswald/MomentOpt.jl/blob/master/src/mopt/Examples.jl
#----------------------------------------------------------------------------------
@testset "Testing smmoptimize with 18d and different magnitude" begin

    # Rermark:
    # It is important NOT to use srand()
    # within the function simulate_empirical_moments!
    # Otherwise, BlackBoxOptim does not find the solution
    #----------------------------------------------------
    # The difference of magniture make it more difficult to find the minimum
    tol18dMean = 0.5

    function functionTest18d(x)

        d = MvNormal(x, eye(18))

        output = OrderedDict{String,Float64}()
        
        draws = rand(d, 1000000)
        for i in 1:18
            output["mean$(i)"] = mean(draws[i,:])
        end

        return output
    end


    t = SMMProblem(options = SMMOptions(maxFuncEvals=10000,saveSteps = 1000))


    #------------------------------------------------------
    dictEmpiricalMoments = OrderedDict{String,Array{Float64,1}}()
    dictEmpiricalMoments["mean1"] = [ -1.0; -1.0]
    dictEmpiricalMoments["mean2"] = [ 0.01; 0.01]
    dictEmpiricalMoments["mean3"] = [ 0.5; 0.5]
    dictEmpiricalMoments["mean4"] = [ -0.5; -0.5]
    dictEmpiricalMoments["mean5"] = [ 0.008; 0.008]
    dictEmpiricalMoments["mean6"] = [-0.7; -0.7]

    dictEmpiricalMoments["mean7"] = [-1.0; -1.0]
    dictEmpiricalMoments["mean8"] = [0.01; 0.01]
    dictEmpiricalMoments["mean9"] = [ 0.5; 0.5]
    dictEmpiricalMoments["mean10"] = [ -0.5; -0.5]
    dictEmpiricalMoments["mean11"] = [ 0.008; 0.008]
    dictEmpiricalMoments["mean12"] = [-0.7; -0.7]

    dictEmpiricalMoments["mean13"] = [-1.0; -1.0]
    dictEmpiricalMoments["mean14"] =  [0.01; 0.01]
    dictEmpiricalMoments["mean15"] = [ 0.5; 0.5]
    dictEmpiricalMoments["mean16"] = [ -0.5; -0.5]
    dictEmpiricalMoments["mean17"] = [ 0.008; 0.008]
    dictEmpiricalMoments["mean18"] = [-0.7; -0.7]

    set_empirical_moments!(t, dictEmpiricalMoments)

    pb = OrderedDict{String,Array{Float64,1}}()

    pb["p1"] = [0.2, 0., 1.]
    pb["p2"] = [0.002, 0., 0.02]
    pb["p3"] = [-0.3,-2., 2.]
    pb["p4"] = [-0.4,-2., 2.]
    pb["p5"] = [0.008,0.,0.02]
    pb["p6"] = [0.4,-2., 2.]

    pb["p7"] = [0.2,0. ,1.]
    pb["p8"] = [0.002,0., 0.02]
    pb["p9"] = [-0.3,-2., 2.]
    pb["p10"] = [-0.4,-2., 2.]
    pb["p11"] = [0.008,0,0.02]
    pb["p12"] = [0.4,-2., 2.]

    pb["p13"] = [0.2, 0., 1.]
    pb["p14"] = [0.002,0,0.02]
    pb["p15"] = [-0.3,-2., 2.]
    pb["p16"] = [-0.4,-2., 2.]
    pb["p17"] = [0.008,0., 0.02]
    pb["p18"] = [0.4,-2., 2.]

    set_priors!(t, pb)

    # A. Set the function: parameter -> simulated moments
    set_simulate_empirical_moments!(t, functionTest18d)

    # B. Construct the objective function, using the function: parameter -> simulated moments
    # and moments' weights:
    construct_objective_function!(t)

    # C. Run the optimization
    # This function first modifies t.bbSetup
    # and then modifies t.bbResults
    smmoptimize!(t, verbose = true)

    @test best_candidate(t.bbResults)[1] ≈ -1.0 atol = tol6dMean
    @test best_candidate(t.bbResults)[2] ≈ 0.01 atol = tol6dMean
    @test best_candidate(t.bbResults)[3] ≈ 0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[4] ≈ -0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[5] ≈ 0.008 atol = tol6dMean
    @test best_candidate(t.bbResults)[6] ≈ -0.7 atol = tol6dMean

    @test best_candidate(t.bbResults)[7] ≈ -1.0 atol = tol6dMean
    @test best_candidate(t.bbResults)[8] ≈ 0.01 atol = tol6dMean
    @test best_candidate(t.bbResults)[9] ≈ 0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[10] ≈ -0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[11] ≈ 0.008 atol = tol6dMean
    @test best_candidate(t.bbResults)[12] ≈ -0.7 atol = tol6dMean

    @test best_candidate(t.bbResults)[13] ≈ -1.0 atol = tol6dMean
    @test best_candidate(t.bbResults)[14] ≈ 0.01 atol = tol6dMean
    @test best_candidate(t.bbResults)[15] ≈ 0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[16] ≈ -0.5 atol = tol6dMean
    @test best_candidate(t.bbResults)[17] ≈ 0.008 atol = tol6dMean
    @test best_candidate(t.bbResults)[18] ≈ -0.7 atol = tol6dMean

end
