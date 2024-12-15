module PUPM
using Ferrite
using Printf
using Tensors

export DynamicParams
include("DynamicParams.jl")

export remove_vtk_files
export get_material_matrix
export assemble_cell!
export assemble_global!
export assemble_external_forces!
export assemble_external_pressure!
export apply_nodal_force!
export calculate_stresses
export calculate_strains
export calculate_strain_energy
export get_material_matrix_derivative_wrt_E
export calculate_cell_volume
export calculate_H
export calculate_average_strain_energy
export LoadCondition
export fem_solver
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
export LoadCondition_3d
export fem_solver_3d
####### optimization functions
export update_upm!
export transfer_to_density!
export transfer_to_young!
#export transform
export filter_density_to_vf!
export upm_update!
export top_upm!
export top_upm_3d!
include("utils.jl")

end
