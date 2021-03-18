# load required codes, and setup the environment -
include("Include.jl")


# -- MAIN ------------------------------------------------------------------------------------------------- #
function sample_bounds_array(path_to_parameter_file::String, objective_function_index_array::Array{Int64,1}, inducer_concentration::Float64)

    # initialize -
    flux_ensemble = zeros(15,1)

    # buld default_parameter_dictionary -
    data_dictionary = generate_parameter_dictionary(path_to_parameter_file)

    # grab the flux bounds array -
    flux_bounds_array = data_dictionary[:flux_bounds_array]
    (number_of_fluxes,number_of_cols) = size(flux_bounds_array)
    flux_bounds_array = (1/1000).*flux_bounds_array
    for flux_index = 1:number_of_fluxes

        # generate the parameter dictionary -
        default_parameter_dictionary = generate_parameter_dictionary(path_to_parameter_file)

        # increase flux index by 1%
        delta = flux_bounds_array[flux_index,2]*0.01
        flux_bounds_array[flux_index,2] = flux_bounds_array[flux_index,2]+delta

        # set the bounds -
        default_parameter_dictionary[:flux_bounds_array] = flux_bounds_array

        # update the bounds array, and the objective function -
        objd = update_objective_function(default_parameter_dictionary,objective_function_index_array)
        objd[:induction_parameter_dictionary][:inducer_concentration] = inducer_concentration
        updated_parameter_dictionary = update_default_bounds_array(objd)

        # ok, get the stuff needed to pass to calculate_optimal_flux_distribution method -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_parameter_dictionary)

        # check the solution quality -
        if (exit_flag == 0 && status_flag == 5)
            flux_ensemble = [flux_ensemble calculated_flux_array]
        end
    end

    return flux_ensemble[:,2:end]
end

function simulate_max_translation(path_to_parameter_file::String, objective_function_index_array::Array{Int64,1},inducer_array::Array{Float64,1})

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
        updated_parameter_dictionary = update_default_bounds_array(objd)

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

# how many simulation points do we have?
number_of_simulation_points = 10

# compute the fluxes -
path_to_parameter_file = "./src/config/Parameters_nominal.json"
objective_index_array = [5] # max translation -

# compute the protein -
inducer_array = collect(exp10.(range(-5,stop=2,length=number_of_simulation_points)))
(dd,fe) = simulate_max_translation(path_to_parameter_file, objective_index_array, inducer_array)

# compute the steady-state protein concentration -
# kdL = dd[:kdL]
# p_array = Float64[]
# translation_rate_array = fe[5,:]
# for translation_rate in translation_rate_array
#     value = (translation_rate/kdL)
#     push!(p_array, value)
# end
