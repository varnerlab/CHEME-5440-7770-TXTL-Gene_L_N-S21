# includes -
include("Include.jl")

function simulate_no_translation_rate(path_to_parameter_file::String, objective_function_index_array::Array{Int64,1},inducer_array::Array{Float64,1})

    # initialize -
    flux_ensemble = zeros(15,1)

    # buld default_parameter_dictionary -
    data_dictionary = generate_parameter_dictionary(path_to_parameter_file)

    # main simulation loop -
    for inducer_concentration in inducer_array

        # generate the parameter dictionary -
        default_parameter_dictionary = generate_parameter_dictionary(path_to_parameter_file)

        # update the bounds array, and the objective function -
        objd = update_objective_function(default_parameter_dictionary,objective_function_index_array)
        objd[:induction_parameter_dictionary][:inducer_concentration] = inducer_concentration
        updated_parameter_dictionary = no_translation_bounds_array(objd)

        # ok, get the stuff needed to pass to calculate_optimal_flux_distribution method -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_parameter_dictionary)

        # check the solution quality -
        if (exit_flag == 0 && status_flag == 5)
            flux_ensemble = [flux_ensemble calculated_flux_array]
        end
    end

    return (data_dictionary, flux_ensemble[:,2:end])
end
# -- MAIN ------------------------------------------------------------------------------------------------- #

# compute the fluxes -
path_to_parameter_file = "./src/config/Parameters_nominal.json"
objective_index_array = [3] # max degradation -

# compute the protein -
inducer_array_no_TL = [100.0]
(dd,fe_no_TL) = simulate_no_translation_rate(path_to_parameter_file, objective_index_array, inducer_array_no_TL)