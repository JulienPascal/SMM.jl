using SMM
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# To make sure Travis CI works as expected
# @test 1 == 2


@testset "testing Types" begin


    @testset "testing SMMOptions" begin

        t = SMMOptions()

        # By default, the optimizer should be :adaptive_de_rand_1_bin_radiuslimited
        #--------------------------------------------------------------------------
        @test t.bbOptimizer == :adaptive_de_rand_1_bin_radiuslimited

    end


    @testset "testing SMMOptions" begin

        t = SMMProblem()

        @test t.iter == 0
        @test typeof(t.priors) == Dict{String, Any} 
        @test typeof(t.empiricalMoments) == Dict{String, Any}
        @test typeof(t.simulatedMoments) == Dict{String, Any}
        @test typeof(t.distanceEmpSimMoments) == Float64
        @test t.simulate_empirical_moments(1.0) == 1.0
        @test t.objective_function(1.0) == 1.0
        @test typeof(t.options) == SMMOptions
        

    end



end