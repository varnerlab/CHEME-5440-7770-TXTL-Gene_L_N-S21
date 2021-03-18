# load required codes, and setup the environment -
include("Include.jl")

function make_binary_stoichiometric_matrix(stoichiometric_matrix::Array{Float64,2})::Array{Float64,2}

    # initialize -
    (number_of_rows, number_of_cols) = size(stoichiometric_matrix)
    binary_stoichiometric_matrix = zeros(number_of_rows, number_of_cols)
    for row_index = 1:number_of_rows
        for col_index = 1:number_of_cols

            # grab value -
            old_value = stoichiometric_matrix[row_index, col_index]
            if (old_value !=0.0)
                binary_stoichiometric_matrix[row_index, col_index] = 1.0
            end
        end
    end

    # return -
    return binary_stoichiometric_matrix
end

function execute_structural_analysis(path_to_parameters_file::String)::NamedTuple

    try

        # buld default_parameter_dictionary -
        data_dictionary = generate_parameter_dictionary(path_to_parameters_file)

        # get the stm -
        stoichiometric_matrix = data_dictionary[:stoichiometric_matrix]

        # compute the binary stm -
        binary_stoichiometric_matrix = make_binary_stoichiometric_matrix(stoichiometric_matrix)

        # compute the metabolite adj array -
        species_adj_array = binary_stoichiometric_matrix*transpose(binary_stoichiometric_matrix)

        # compute the reaction adj array -
        reaction_adj_array = transpose(binary_stoichiometric_matrix)*binary_stoichiometric_matrix

        # setup the named return tuple -
        return_tuple = (A_m=species_adj_array, A_r=reaction_adj_array, stoichiometric_matrix=stoichiometric_matrix, binary_stoichiometric_matrix=binary_stoichiometric_matrix)

        # return -
        return return_tuple

    catch error
        throw(error)
    end
end


# load the stm matrix -
path_to_parameter_file = "./src/config/Parameters_nominal.json"

# do some structurual analysis -
structural_analysis_tuple = execute_structural_analysis(path_to_parameter_file)