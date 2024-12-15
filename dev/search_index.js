var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Pages = [\"functions.md\"]","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"Pages = [\"functions.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = PUPM","category":"page"},{"location":"#PUPM","page":"Home","title":"PUPM","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"welcome to my Documentation for PUPM.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PUPM]","category":"page"},{"location":"#PUPM.apply_nodal_force!-NTuple{4, Any}","page":"Home","title":"PUPM.apply_nodal_force!","text":"Function to apply nodal forces to external force vector\n\napply_nodal_force!(grid, node_set_name, load_vector, f_ext)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.apply_nodal_force_3d!-NTuple{4, Any}","page":"Home","title":"PUPM.apply_nodal_force_3d!","text":"function to apply nodal forces to external force vector in 3D\n\napply_nodal_force_3d!(grid, node_set_name, load_vector, f_ext)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_cell!-NTuple{4, Any}","page":"Home","title":"PUPM.assemble_cell!","text":"Function for the local stiffness matrix(element stiffness matrix)\n\nassemble_cell!(ke, cell_values, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_cell_3d!-NTuple{4, Any}","page":"Home","title":"PUPM.assemble_cell_3d!","text":"function to assemble the local stiffness matrix for 3D\n\nassemble_cell_3d!(ke, cell_values, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_external_forces!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_external_forces!","text":"Function for  external forces from surface tractions\n\nassemble_external_forces!(f_ext, dh, facetset, facet_values, traction)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_external_forces_3d!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_external_forces_3d!","text":"function to assemble external forces from surface tractions in 3D\n\nassemble_external_forces_3d!(f_ext, dh, facetset, facet_values, traction)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_external_pressure!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_external_pressure!","text":"Function for  external forces from pressure\n\nassemble_external_pressure!(f_ext, dh, facetset, facet_values, pressure)\n\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_external_pressure_3d!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_external_pressure_3d!","text":"function to assemble external forces from pressure in 3D\n\nassemble_external_pressure_3d!(f_ext, dh, facetset, facet_values, pressure)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_global!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_global!","text":"Function for the global stiffness matrix\n\nassemble_global!(K, dh, cell_values, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.assemble_global_3d!-NTuple{5, Any}","page":"Home","title":"PUPM.assemble_global_3d!","text":"function to assemble the global stiffness matrix for 3D\n\nassemble_global_3d!(K, dh, cell_values, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_H-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_H","text":"Function to calculate parameter H\n\ncalculate_H(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_H_3d-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_H_3d","text":"function to calculate parameter H in 3D\n\ncalculate_H_3d(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_average_strain_energy-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_average_strain_energy","text":"Function to calculate the average strain energy\n\ncalculate_average_strain_energy(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_cell_volume-Tuple{Any}","page":"Home","title":"PUPM.calculate_cell_volume","text":"Function to calculate the cell volume\n\ncalculate_cell_volume(cv)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_strain_energy-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_strain_energy","text":"Function to calculate element strain energy\n\ncalculate_strain_energy(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_strain_energy_3d-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_strain_energy_3d","text":"function to calculate element strain energy in 3D\n\ncalculate_strain_energy_3d(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_strains-NTuple{4, Any}","page":"Home","title":"PUPM.calculate_strains","text":"Function to calculate strains\n\ncalculate_strains(grid, dh, cv, u)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_strains_3d-NTuple{4, Any}","page":"Home","title":"PUPM.calculate_strains_3d","text":"function to calculate strains in 3D\n\ncalculate_strains_3d(grid, dh, cv, u)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_stresses-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_stresses","text":"Function to calculate stresses\n\ncalculate_stresses(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.calculate_stresses_3d-NTuple{6, Any}","page":"Home","title":"PUPM.calculate_stresses_3d","text":"function to calculate stresses in 3D\n\ncalculate_stresses_3d(grid, dh, cv, u, E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.fem_solver-Tuple{DynamicParams}","page":"Home","title":"PUPM.fem_solver","text":"Main finite element solver for two dimension\n\nfem_solver(par::DynamicParams)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.fem_solver_3d-Tuple{DynamicParams}","page":"Home","title":"PUPM.fem_solver_3d","text":"Main finite element solver for 3D\n\nfem_solver_3d(par::DynamicParams)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.filter_density_to_vf!-NTuple{6, Any}","page":"Home","title":"PUPM.filter_density_to_vf!","text":"function to apply volume fraction\n\nexample:\n\nρnew = rand(Float64, 8)\nnx, ny , nz = 2 , 2 , 2\nvf = 0.5\nη = π/4\n\nρ =  filter_density_to_vf!(ρnew, vf, nx, ny, nz, η)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.get_material_matrix-Tuple{Any, Any}","page":"Home","title":"PUPM.get_material_matrix","text":"Get the material matrix C\n\nget_material_matrix(E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.get_material_matrix_3d-Tuple{Any, Any}","page":"Home","title":"PUPM.get_material_matrix_3d","text":"function to get the material matrix for 3D\n\nget_material_matrix_3d(E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.get_material_matrix_derivative_wrt_E-Tuple{Any, Any}","page":"Home","title":"PUPM.get_material_matrix_derivative_wrt_E","text":"Function to calculate the derivative of the strain energy with respect to the Young's modulus\n\nget_material_matrix_derivative_wrt_E(E, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.get_material_matrix_derivative_wrt_E_3d-Tuple{Any, Any}","page":"Home","title":"PUPM.get_material_matrix_derivative_wrt_E_3d","text":"function to get the derivative of the material matrix with respect to the Young's modulus in 3D\n\nget_material_matrix_derivative_wrt_E_3d(E::Vector{Float64,1}, ν)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.remove_vtk_files-Tuple{String}","page":"Home","title":"PUPM.remove_vtk_files","text":"function to remove vtu file \n\nremove_vtk_files(directory::String)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.top_upm!-Tuple{DynamicParams, String, String}","page":"Home","title":"PUPM.top_upm!","text":"function to perform topology optimization using UPM approach (2D case)\n\ntop_upm!(par::DynamicParams, name_of_file::String, directory::String)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.top_upm_3d!-Tuple{DynamicParams, String, String}","page":"Home","title":"PUPM.top_upm_3d!","text":"function to perform topology optimization using UPM approach (3D case)\n\ntop_upm_3d!(par::DynamicParams, name_of_file::String, directory::String)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.transfer_to_density!-Tuple{Vector{Float64}, Float64, Float64, Int64}","page":"Home","title":"PUPM.transfer_to_density!","text":"function to convert Young's modulus to density     E = E0(rho/rho0)^gamma\n\nExample:\n\nEnew = rand(Float64, 5)\nE0 = 1.0\nρ0 = 1.0\nγ = 1\nρ = transfer_to_density!(Enew, E0, ρ0, γ)\n\n\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.transfer_to_young!-Tuple{Vector{Float64}, Float64, Float64, Int64, Float64, Float64}","page":"Home","title":"PUPM.transfer_to_young!","text":"function to transfer density to young modulus example: ρnew = rand(Float64, 4) E0 = 1.0 ρ0 = 1.0 γ = 1 Emin = 1e-4 Emax = 1.0\n\ntransfertoyoung!(ρnew , E0, ρ0, γ, Emin , Emax)\n\n\n\n\n\n","category":"method"},{"location":"#PUPM.update_upm!-Tuple{Int64, Vector{Float64}, Vector{Float64}, Float64, Float64}","page":"Home","title":"PUPM.update_upm!","text":"function for update Young's modulus based on UPM approach\n\nExample:\n\nk = 1\nE = rand(Float64, 5)\nH = rand(Float64, 5)\nEmax = 1.0\nEmin = 1e-4\n# test code\nEnew = update_upm!(k, E, H,Emax,Emin)\n\n\n\n\n\n\n","category":"method"}]
}
