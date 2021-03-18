# Load some packages -
using CodeGenerator
using DelimitedFiles

# What files am I going to load/write?
path_to_reaction_file = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/classes/CHEME-7770-S19/code/txtl_prelim_model/config/Reactions.net"
path_to_debug_file = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/classes/CHEME-7770-S19/code/txtl_prelim_model/config/Debug.txt"
path_to_stochiometric_matrix = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/classes/CHEME-7770-S19/code/txtl_prelim_model/config/Network.dat"
path_to_bounds_matrix = "/Users/jeffreyvarner/Desktop/Desktop_Jeffrey_MacBook_Pro/classes/CHEME-7770-S19/code/txtl_prelim_model/config/Bounds.dat"

# Build the st matrix -
stoichiometrix_matrix = build_stoichiometric_matrix_from_vff_file(path_to_reaction_file)
writedlm(path_to_stochiometric_matrix, stoichiometrix_matrix)

# Build the bounds matrix -
flux_bounds_matrix = build_reaction_bounds_matrix_from_vff_file(path_to_reaction_file)
writedlm(path_to_bounds_matrix, flux_bounds_matrix)

# build the debug file to check generated code -
write_debug_report(path_to_reaction_file, path_to_debug_file)
