module PUPM
using Ferrite
using Printf
using Tensors
using Statistics
using LinearAlgebra
using ForwardDiff

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
export fem_solver_combine
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
export fem_solver_3d_combine
################################
include("function_fem_3d.jl")
export update_upm
export transfer_to_density
export transfer_to_young
export filter_density_to_vf!
#export upm_update
export top_upm
export top_upm_3d

export top_upm_combine
export top_upm_3d_combine
export optim_2D_combine
export optim_3D_combine
include("optim_function.jl")
################################
####### optimization functions
export remove_vtk_files

include("utils.jl")
################################
export C_orthotropic
export compute_derivatives
export assemble_cell_orthotropic!
export assemble_global_orthotropic!
export calculate_strain_energy_orthotropic
export dC_orthotropic_dEx
export dC_orthotropic_dEy
export calculate_Hx
export calculate_Hy
export FEMSolver_orthotropic
export fem_solver_combine_orthotropic
export compute_gamma
export transfer_to_young_y_dir
export top_upm_orthotropic
export optim_2D_combine_orthotropic
include("functions_orthopedic jl")
################################


end
