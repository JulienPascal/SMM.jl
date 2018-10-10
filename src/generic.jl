"""
  set_simulate_empirical_moments!(SMMProblem, f::Function)

Function to set the field simulate_empirical_moments for an SMMProblem.
The function simulate_empirical_moments takes parameter values and return
the corresponding simulate moments values.

"""
function set_simulate_empirical_moments!(SMMProblem, f::Function)

  SMMProblem.simulate_empirical_moments = f

end

"""
  construct_objective_function!(SMMProblem, f::Function)

Function that construct an objective function, using the function
SMMProblem.simulate_empirical_moments. For the moment, the objective function
returns the mean of the percentage deviation of simulated moments from
their true value.

"""
function construct_objective_function!(SMMProblem, objectiveType::Symbol = :percent)

  if objectiveType == :percent

    function objective_function(params::Vector)

      simulatedMoments, convegence = try

          SMMProblem.simulate_empirical_moments(params), 1

      catch errorSimulation

        info("$(errorSimulation) with parameter values = $(params)")

          Dict{String, Any}(), 0

      end

      return simulatedMoments, convegence

    end

  else

    error("Error in construct_objective_function. objectiveType currently only support :percent.")

  end


  SMMProblem.objective_function = objective_function


end
