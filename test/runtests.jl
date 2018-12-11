# Let's increase the number of workers
#-------------------------------------
maxNbWorkers = 3
while nworkers() < maxNbWorkers
  addprocs(1)
end

@everywhere using SMM
@everywhere using Distributions
using Optim

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@static if VERSION < v"0.7.0"
    @everywhere using DataStructures
else
    @everywhere using OrderedCollections
    # Because of this issue (https://github.com/JuliaIO/JLD2.jl/issues/107)
    # we also need to import BlackBoxOptim and Optim to load and save
    #----------------------------------------------------------------------
    # (no need when using Julia 0.6)
    using BlackBoxOptim

end


# To make sure Travis CI works as expected
# @test 1 == 2


@testset "testing Types" begin


    @testset "testing SMMOptions" begin

        t = SMMOptions()

        # By default, the optimizer should be :adaptive_de_rand_1_bin_radiuslimited
        #--------------------------------------------------------------------------
        @test t.globalOptimizer == :dxnes
	      @test t.localOptimizer == :LBFGS

    end


    @testset "testing SMMOptions" begin

        t = SMMProblem()

        # When initialized, iter is equal to 0
        @test t.iter == 0
        @test typeof(t.priors) == OrderedDict{String,Array{Float64,1}}
        @test typeof(t.empiricalMoments) == OrderedDict{String,Array{Float64,1}}
        @test typeof(t.simulatedMoments) == OrderedDict{String, Float64}
        @test typeof(t.distanceEmpSimMoments) == Float64
        # the functions t.simulate_empirical_moments are initialized with x->x
        @test t.simulate_empirical_moments(1.0) == 1.0
        @test t.objective_function(1.0) == 1.0
        @test typeof(t.options) == SMMOptions


    end



end


@testset "testing creation of Cartesian grids" begin

  a = zeros(3)
  b = ones(3)
  nums = 2

  points = cartesian_grid(a, b, nums)

  @test points[1] == (0.0, 0.0, 0.0)
  @test points[2] == (1.0, 0.0, 0.0)
  @test points[3] == (0.0, 1.0, 0.0)
  @test points[4] == (1.0, 1.0, 0.0)
  @test points[5] == (0.0, 0.0, 1.0)
  @test points[6] == (1.0, 0.0, 1.0)
  @test points[7] == (0.0, 1.0, 1.0)
  @test points[8] == (1.0, 1.0, 1.0)

  points = create_grid(a, b, nums)

  @test points[1,:] == [0.0; 0.0; 0.0]
  @test points[2,:] == [1.0, 0.0, 0.0]
  @test points[3,:] == [0.0, 1.0, 0.0]
  @test points[4,:] == [1.0, 1.0, 0.0]
  @test points[5,:] == [0.0, 0.0, 1.0]
  @test points[6,:] == [1.0, 0.0, 1.0]
  @test points[7,:] == [0.0, 1.0, 1.0]
  @test points[8,:] == [1.0, 1.0, 1.0]


end

@testset "testing creation of random grids" begin

  a = zeros(3)
  b = [1.0; 2.0; 3.0]
  numPoints = 100


  @testset "drawing from a multivariate normal" begin

    points = create_grid_stochastic(a, b, numPoints, gridType = :normal, alpha = 0.05)

    @test size(points, 1) == numPoints
    @test size(points, 2) == length(a)

    for i=1:size(points, 1)
      for j=1:size(points, 2)
        @test points[i,j] >= a[j]
        @test points[i,j] <= b[j]
      end
    end

  end

  @testset "drawing from a multivariate uniform" begin

    points = create_grid_stochastic(a, b, numPoints, gridType = :uniform)

    @test size(points, 1) == numPoints
    @test size(points, 2) == length(a)

    for i=1:size(points, 1)
      for j=1:size(points, 2)
        @test points[i,j] >= a[j]
        @test points[i,j] <= b[j]
      end
    end

  end


end

@testset "testing checks on global algo" begin

    listValidGlobalOptimizers = [:dxnes, :adaptive_de_rand_1_bin_radiuslimited, :xnes,
                             :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin,
                             :generating_set_search, :de_rand_1_bin,
                             :separable_nes, :resampling_inheritance_memetic_search,
                             :probabilistic_descent, :resampling_memetic_search,
                             :de_rand_2_bin_radiuslimited, :de_rand_2_bin,
                             :random_search, :simultaneous_perturbation_stochastic_approximation]


    for globalOptim in listValidGlobalOptimizers

      @test is_global_optimizer(globalOptim) == true

    end

end


@testset "testing checks on local algo" begin

  listValidLocalOptimizers = [:NelderMead, :SimulatedAnnealing, :ParticleSwarm,
                            :BFGS, :LBFGS]


    for globalOptim in listValidLocalOptimizers

      @test is_local_optimizer(globalOptim) == true

    end

end


@testset "testing loading priors and empirical moments" begin


   @testset "testing read_priors" begin

        dictPriors = read_priors(joinpath(Pkg.dir("SMM"), "test/priorsTest.csv"))

        @test typeof(dictPriors) == OrderedDict{String,Array{Float64,1}}
        # First component stores the value
        @test dictPriors["alpha"][1] == 0.5
        # Second component stores the lower bound:
        @test dictPriors["alpha"][2] == 0.01
        # Third component stores the upper bound:
        @test dictPriors["alpha"][3] == 0.9


   end

   @testset "testing read_empirical_moments" begin

        dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))

        @test typeof(dictEmpiricalMoments) == OrderedDict{String,Array{Float64,1}}

        # First component stores the value
        @test dictEmpiricalMoments["meanU"][1] == 0.05
        # Second component stores the weight associated to
        @test dictEmpiricalMoments["meanU"][2] == 0.05


   end


end


@testset "set_priors!, set_empirical_moments!" begin

    t = SMMProblem()

    @testset "testing set_priors!" begin

        dictPriors = read_priors(joinpath(Pkg.dir("SMM"), "test/priorsTest.csv"))

        set_priors!(t, dictPriors)

        @test t.priors == dictPriors

    end

    @testset "set_empirical_moments!" begin

        dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))

        set_empirical_moments!(t, dictEmpiricalMoments)

        @test t.empiricalMoments == dictEmpiricalMoments

    end
end


@testset "testing the construction of the objective function" begin


   @testset "set_simulate_empirical_moments!" begin

        function functionTest(x::Vector)

            output = OrderedDict{String,Float64}()
            output["mom1"] = x[1]
            output["mom2"] = x[2]

            return output
        end

        t = SMMProblem()

        set_simulate_empirical_moments!(t, functionTest)

        x1Value = 1.0
        x2Value = 2.0
        simulatedMoments = t.simulate_empirical_moments([x1Value; x2Value])
        @test simulatedMoments["mom1"]  == x1Value
        @test simulatedMoments["mom2"]  == x2Value

   end

   @testset "Testing construct_objective_function!" begin

        srand(1234)
        tol1dMean = 0.01

        function functionTest(x::Vector)

            output = OrderedDict{String,Float64}()
            d = Normal(x[1])
            output["meanU"] = mean(rand(d, 100000))

            return output
        end

        t = SMMProblem()

        set_simulate_empirical_moments!(t, functionTest)

        # For the test to make sense, we need to set the field
        # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
        #------------------------------------------------------
        dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))
        set_empirical_moments!(t, dictEmpiricalMoments)

        # A. Set the function: parameter -> simulated moments
        set_simulate_empirical_moments!(t, functionTest)

        # B. Construct the objective function, using the function: parameter -> simulated moments
        # and moments' weights:
        construct_objective_function!(t)

        # The objective function should be very close to 0 when
        # evaluated at the true value (modulo randomness)
        @test t.objective_function([dictEmpiricalMoments["meanU"][1]]) ≈ 0. atol = tol1dMean


    end

    @testset "Testing generate_bbSearchRange" begin

        # A.
        #----
        t = SMMProblem()

        dictPriors = read_priors(joinpath(Pkg.dir("SMM"), "test/priorsTest.csv"))

        set_priors!(t, dictPriors)

        testSearchRange = generate_bbSearchRange(t)

        @test testSearchRange[1][1] == 0.01
        @test testSearchRange[1][2] == 0.9
        @test testSearchRange[2][1] == 0.0
        @test testSearchRange[2][2] == 1.0

        # B.
        #---
        t = SMMProblem()

        dictPriors = OrderedDict{String,Array{Float64,1}}()
        dictPriors["mu1"] = [0., -5.0, 5.0]
        dictPriors["mu2"] = [0., -15.0, -10.0]
        dictPriors["mu3"] = [0., -20.0, -15.0]

        set_priors!(t, dictPriors)

        testSearchRange = generate_bbSearchRange(t)

        @test testSearchRange[1][1] == -5.0
        @test testSearchRange[1][2] == 5.0
        @test testSearchRange[2][1] == -15.0
        @test testSearchRange[2][2] == -10.0
        @test testSearchRange[3][1] == -20.0
        @test testSearchRange[3][2] == -15.0


    end


end

@testset "testing loading and saving an optimization" begin


    srand(1234)
    tol1dMean = 0.01

    function functionTest(x::Vector)

        output = OrderedDict{String,Float64}()
        d = Normal(x[1])
        output["meanU"] = mean(rand(d, 100000))

        return output
    end

    t = SMMProblem()

    set_simulate_empirical_moments!(t, functionTest)

    # For the test to make sense, we need to set the field
    # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
    #------------------------------------------------------
    dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))
    set_empirical_moments!(t, dictEmpiricalMoments)

    # A. Set the function: parameter -> simulated moments
    set_simulate_empirical_moments!(t, functionTest)

    # B. Construct the objective function, using the function: parameter -> simulated moments
    # and moments' weights:
    construct_objective_function!(t)


    saveSMMOptim(t, saveName = "iamatest")
    t2 = loadSMMOptim("iamatest")

    # Test the objective function is correctly loaded
    #------------------------------------------------
    # [TODO] test other fields
    @test t.objective_function([dictEmpiricalMoments["meanU"][1]]) ≈ t2.objective_function([dictEmpiricalMoments["meanU"][1]]) atol = tol1dMean

end


@testset "Testing smmoptimize" begin


    # 1d problem
    #-----------
    @testset "Testing smmoptimize with 1d" begin

        # Rermark:
        # It is important NOT to use srand()
        # within the function simulate_empirical_moments!
        # Otherwise, BlackBoxOptim does not find the solution
        #----------------------------------------------------
        tol1dMean = 0.1

        @everywhere function functionTest1d(x)

            d = Normal(x[1])
            output = OrderedDict{String,Float64}()

            output["meanU"] = mean(rand(d, 1000000))

            return output
        end


        t = SMMProblem(options = SMMOptions(maxFuncEvals=200,saveSteps = 100))

        # For the test to make sense, we need to set the field
        # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
        #------------------------------------------------------
        dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))
        set_empirical_moments!(t, dictEmpiricalMoments)


        dictPriors = OrderedDict{String,Array{Float64,1}}()
        dictPriors["mu1"] = [0., -2.0, 2.0]
        set_priors!(t, dictPriors)

        # A. Set the function: parameter -> simulated moments
        #----------------------------------------------------
        set_simulate_empirical_moments!(t, functionTest1d)

        # B. Construct the objective function, using the function: parameter -> simulated moments
        # and moments' weights:
        #----------------------------------------------------
        construct_objective_function!(t)

        smmoptimize!(t, verbose = true)

        @test best_candidate(t.bbResults)[1] ≈ 0.05 atol = tol1dMean

        # C. Testing refinement of the global max using a local routine
        #--------------------------------------------------------------
        @test smm_refine_globalmin!(t, verbose = true)[1] ≈ 0.05 atol = tol1dMean

        @test smm_local_minimizer(t)[1] ≈ 0.05 atol = tol1dMean

    end

    @testset "Testing local to global 1d" begin


      tol1dMean = 0.1

      @everywhere function functionTest1d(x)

          # When using one of the deterministic methods of Optim,
          # we can safely "control" for randomness
          #-----------------------------------------------------
          srand(1234)
          d = Normal(x[1])
          output = OrderedDict{String,Float64}()

          output["meanU"] = mean(rand(d, 1000000))

          return output
      end


      t = SMMProblem(options = SMMOptions(maxFuncEvals=200,saveSteps = 100))

      # For the test to make sense, we need to set the field
      # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
      #------------------------------------------------------
      dictEmpiricalMoments = read_empirical_moments(joinpath(Pkg.dir("SMM"), "test/empiricalMomentsTest.csv"))
      set_empirical_moments!(t, dictEmpiricalMoments)


      dictPriors = OrderedDict{String,Array{Float64,1}}()
      dictPriors["mu1"] = [0., -2.0, 2.0]
      set_priors!(t, dictPriors)

      # A. Set the function: parameter -> simulated moments
      #----------------------------------------------------
      set_simulate_empirical_moments!(t, functionTest1d)

      # B. Construct the objective function, using the function: parameter -> simulated moments
      # and moments' weights:
      #----------------------------------------------------
      construct_objective_function!(t)

      local_to_global!(t, nums = maxNbWorkers, verbose = true)

      @test smm_local_minimum(t) ≈ 0.0 atol = tol1dMean

      @test smm_local_minimizer(t)[1] ≈ 0.05 atol = tol1dMean

    end


    # 2d problem
    #-----------
    @testset "Testing smmoptimize with 2d and same magnitude" begin

        # Rermark:
        # It is important NOT to use srand()
        # within the function simulate_empirical_moments!
        # Otherwise, BlackBoxOptim does not find the solution
        #----------------------------------------------------
        tol2dMean = 0.2

        @everywhere function functionTest2d(x)

            d = MvNormal( [x[1]; x[2]], eye(2))
            output = OrderedDict{String,Float64}()

            draws = rand(d, 1000000)
            output["mean1"] = mean(draws[1,:])
            output["mean2"] = mean(draws[2,:])

            return output
        end


        t = SMMProblem(options = SMMOptions(maxFuncEvals=1000,saveSteps = 500))


        # For the test to make sense, we need to set the field
        # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
        #------------------------------------------------------
        dictEmpiricalMoments = OrderedDict{String,Array{Float64,1}}()
        dictEmpiricalMoments["mean1"] = [1.0; 1.0]
        dictEmpiricalMoments["mean2"] = [-1.0; -1.0]
        set_empirical_moments!(t, dictEmpiricalMoments)


        dictPriors = OrderedDict{String,Array{Float64,1}}()
        dictPriors["mu1"] = [0., -5.0, 5.0]
        dictPriors["mu2"] = [0., -5.0, 5.0]
        set_priors!(t, dictPriors)

        # A. Set the function: parameter -> simulated moments
        set_simulate_empirical_moments!(t, functionTest2d)

        # B. Construct the objective function, using the function: parameter -> simulated moments
        # and moments' weights:
        construct_objective_function!(t)

        # C. Run the optimization
        # This function first modifies t.bbSetup
        # and then modifies t.bbResults
        smmoptimize!(t, verbose = true)

        @test best_candidate(t.bbResults)[1] ≈ 1.0 atol = tol2dMean
        @test best_candidate(t.bbResults)[2] ≈ -1.0 atol = tol2dMean



    end

    # 2d problem
    #-----------
    @testset "Testing local_to_global! with 2d and same magnitude" begin

        # Rermark:
        # It is important NOT to use srand()
        # within the function simulate_empirical_moments!
        # Otherwise, BlackBoxOptim does not find the solution
        #----------------------------------------------------
        tol2dMean = 0.2

        @everywhere function functionTest2d(x)

            d = MvNormal( [x[1]; x[2]], eye(2))
            output = OrderedDict{String,Float64}()

            # When using one of the deterministic methods of Optim,
            # we can safely "control" for randomness
            #-----------------------------------------------------
            srand(1234)
            draws = rand(d, 1000000)
            output["mean1"] = mean(draws[1,:])
            output["mean2"] = mean(draws[2,:])

            return output
        end


        t = SMMProblem(options = SMMOptions(maxFuncEvals=1000,saveSteps = 500))


        # For the test to make sense, we need to set the field
        # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
        #------------------------------------------------------
        dictEmpiricalMoments = OrderedDict{String,Array{Float64,1}}()
        dictEmpiricalMoments["mean1"] = [1.0; 1.0]
        dictEmpiricalMoments["mean2"] = [-1.0; -1.0]
        set_empirical_moments!(t, dictEmpiricalMoments)


        dictPriors = OrderedDict{String,Array{Float64,1}}()
        dictPriors["mu1"] = [0., -5.0, 5.0]
        dictPriors["mu2"] = [0., -5.0, 5.0]
        set_priors!(t, dictPriors)

        # A. Set the function: parameter -> simulated moments
        set_simulate_empirical_moments!(t, functionTest2d)

        # B. Construct the objective function, using the function: parameter -> simulated moments
        # and moments' weights:
        construct_objective_function!(t)

        local_to_global!(t, nums = nworkers(), verbose = true)

        @test smm_local_minimum(t) ≈ 0.0 atol = tol2dMean

        @test smm_local_minimizer(t)[1] ≈ 1.0 atol = tol2dMean
        @test smm_local_minimizer(t)[2] ≈ - 1.0 atol = tol2dMean

      end

    # 2d problem
    #-----------
    @testset "Testing smmoptimize with 2d with a 1-order magnitude difference" begin

        # Rermark:
        # It is important NOT to use srand()
        # within the function simulate_empirical_moments!
        # Otherwise, BlackBoxOptim does not find the solution
        #----------------------------------------------------
        # The difference of magniture make it more difficult to find the minimum
        tol2dMean = 0.5

        function functionTest2d(x)

            d = MvNormal( [x[1]; x[2]], eye(2))
            output = OrderedDict{String,Float64}()

            draws = rand(d, 1000000)
            output["mean1"] = mean(draws[1,:])
            output["mean2"] = mean(draws[2,:])

            return output
        end


        t = SMMProblem(options = SMMOptions(maxFuncEvals=1000,saveSteps = 500))


        # For the test to make sense, we need to set the field
        # t.empiricalMoments::OrderedDict{String,Array{Float64,1}}
        #------------------------------------------------------
        dictEmpiricalMoments = OrderedDict{String,Array{Float64,1}}()
        dictEmpiricalMoments["mean1"] = [ 1.0; 1.0]
        dictEmpiricalMoments["mean2"] = [-12.0; 12.0]
        set_empirical_moments!(t, dictEmpiricalMoments)


        dictPriors = OrderedDict{String,Array{Float64,1}}()
        dictPriors["mu1"] = [0., -5.0, 5.0]
        dictPriors["mu2"] = [0., -15.0, -10.0]
        set_priors!(t, dictPriors)

        # A. Set the function: parameter -> simulated moments
        set_simulate_empirical_moments!(t, functionTest2d)

        # B. Construct the objective function, using the function: parameter -> simulated moments
        # and moments' weights:
        construct_objective_function!(t)

        # C. Run the optimization
        # This function first modifies t.bbSetup
        # and then modifies t.bbResults
        smmoptimize!(t, verbose = true)

        @test best_candidate(t.bbResults)[1] ≈  1.0 atol = tol2dMean
        @test best_candidate(t.bbResults)[2] ≈ -12.0 atol = tol2dMean

    end


end
