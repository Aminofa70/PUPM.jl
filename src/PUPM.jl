module PUPM
using Ferrite
using Printf
using Tensors
using Statistics


export DynamicParams
include("DynamicParams.jl")

export get_material_matrix
export assemble_cell!
export assemble_global!
export assemble_external_forces!
export assemble_external_pressure!
export vertexdofs
export nodeid_to_vertexindex
export apply_nodal_force!
export calculate_stresses
export calculate_strains
export calculate_strain_energy
export get_material_matrix_derivative_wrt_E
export calculate_cell_volume
export calculate_H
export calculate_average_strain_energy
export compute_nodal_data
export LoadCondition
export fem_solver
include("function_fem_2d.jl")

################################
#### three dimensional functions for fem
export apply_nodal_force_3d!
export get_material_matrix_3d
export assemble_cell_3d!
export assemble_global_3d!
export assemble_external_forces_3d!
export assemble_external_pressure_3d!
export calculate_stresses_3d
export calculate_strains_3d
export calculate_strain_energy_3d
export get_material_matrix_derivative_wrt_E_3d
export calculate_cell_volume_3d
export calculate_H_3d
export compute_nodal_data_3D
export LoadCondition_3d
export fem_solver_3d
################################
include("function_fem_3d.jl")
export update_upm
export transfer_to_density
export transfer_to_young
export filter_density_to_vf
export upm_update
export top_upm
export top_upm_3d
include("optim_function.jl")
################################
####### optimization functions
export remove_vtk_files

include("utils.jl")

end
