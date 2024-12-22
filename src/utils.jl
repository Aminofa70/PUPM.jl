"""
function to calculate the mean of an array of numbers
    ```
    calculate_mean(numbers::Vector{<:Number}) 
    ```

"""
function calculate_mean(numbers::Vector{<:Number})
    # Check if the array is empty to avoid errors
    if isempty(numbers)
        error("Cannot calculate mean of an empty array")
    end
    # Calculate the mean
    sum(numbers) / length(numbers)
end

"""


"""
function standard_deviation(data::AbstractVector{<:Number})
    # Ensure the data is not empty
    @assert length(data) > 1 "Data vector must have at least two elements."
    
    # Compute the mean
    m = calculate_mean(data)
    
    # Compute the sum of squared differences from the mean
    sum_sq_diff = sum((x - m)^2 for x in data)
    
    # Compute the sample standard deviation (note the division by n-1)
    return sqrt(sum_sq_diff / (length(data) - 1))
end

"""
function to remove vtu file 
```
remove_vtk_files(directory::String)
```
"""
function remove_vtk_files(directory::String)
    files = readdir(directory)
    for file in files
        if endswith(file, ".vtu")
            filepath = joinpath(directory, file)
            rm(filepath, force = true)
        end
    end
end
############################################
"""
Get the material matrix C
```
get_material_matrix(E, ν)
```
"""
function get_material_matrix(E, ν)
    C_voigt = E * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end

# Function to assemble the local stiffness matrix
"""
Function for the local stiffness matrix(element stiffness matrix)
```
assemble_cell!(ke, cell_values, E, ν)
```
"""
function assemble_cell!(ke, cell_values, E, ν)
    C = get_material_matrix(E, ν)
    for qp in 1:getnquadpoints(cell_values)
        dΩ = getdetJdV(cell_values, qp)
        for i in 1:getnbasefunctions(cell_values)
            ∇Ni = shape_gradient(cell_values, qp, i)
            for j in 1:getnbasefunctions(cell_values)
                ∇δNj = shape_symmetric_gradient(cell_values, qp, j)
                ke[i, j] += (∇Ni ⊡ C ⊡ ∇δNj) * dΩ
            end
        end
    end
    return ke
end
# 
"""
Function for the global stiffness matrix
```
assemble_global!(K, dh, cell_values, E, ν)
```
"""
function assemble_global!(K, dh, cell_values, E, ν)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)

        # Check if E is scalar or vector
        local_E = isa(E, AbstractVector) ? E[cell_index] : E

        assemble_cell!(ke, cell_values, local_E, ν)
        assemble!(assembler, celldofs(cell), ke)
    end

    return K
end

# function assemble_global!(K, dh, cell_values, E, ν)
#     n_basefuncs = getnbasefunctions(cell_values)
#     ke = zeros(n_basefuncs, n_basefuncs)
#     assembler = start_assemble(K)
#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cell_values, cell)
#         fill!(ke, 0.0)
#         assemble_cell!(ke, cell_values, E[cell_index], ν)
#         assemble!(assembler, celldofs(cell), ke)
#     end
#     return K
# end
# Function to assemble external forces from surface tractions
"""
Function for  external forces from surface tractions
```
assemble_external_forces!(f_ext, dh, facetset, facet_values, traction)
```
"""
function assemble_external_forces!(f_ext, dh, facetset, facet_values, traction)
    fe_ext = zeros(getnbasefunctions(facet_values))
    for facet in FacetIterator(dh, facetset)
        reinit!(facet_values, facet)
        fill!(fe_ext, 0.0)
        for qp in 1:getnquadpoints(facet_values)
            dΓ = getdetJdV(facet_values, qp)
            for i in 1:getnbasefunctions(facet_values)
                fe_ext[i] += traction ⋅ shape_value(facet_values, qp, i) * dΓ
            end
        end
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
# Function to assemble external forces from pressure
"""
Function for  external forces from pressure
```
assemble_external_pressure!(f_ext, dh, facetset, facet_values, pressure)

```
"""
function assemble_external_pressure!(f_ext, dh, facetset, facet_values, pressure)
    fe_ext = zeros(getnbasefunctions(facet_values))
    for facet in FacetIterator(dh, facetset)
        reinit!(facet_values, facet)
        fill!(fe_ext, 0.0)
        for qp in 1:getnquadpoints(facet_values)
            p = pressure * getnormal(facet_values, qp)
            dΓ = getdetJdV(facet_values, qp)
            for i in 1:getnbasefunctions(facet_values)
                fe_ext[i] += p ⋅ shape_value(facet_values, qp, i) * dΓ
            end
        end
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
# Function to apply nodal forces to external force vector
"""
Function to apply nodal forces to external force vector
```
apply_nodal_force!(grid, node_set_name, load_vector, f_ext)
```
"""
function apply_nodal_force!(grid, node_set_name, load_vector, f_ext)
    node_set = Ferrite.getnodeset(grid, node_set_name)
    for node_id in node_set
        f_ext[2*node_id-1] += load_vector[1]  # x-component
        f_ext[2*node_id] += load_vector[2]      # y-component
    end
    return f_ext
end
# Function to calculate stresses
"""
Function to calculate stresses
```
calculate_stresses(grid, dh, cv, u, E, ν)
```
"""
function calculate_stresses(grid, dh, cv, u, E, ν)
    qp_stresses = [
        [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
        for _ in 1:getncells(grid)]
    avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:3)...)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Check if E is scalar or vector
        local_E = isa(E, AbstractVector) ? E[cell_index] : E

        C = get_material_matrix(local_E, ν)
        cell_stresses = qp_stresses[cellid(cell)]

        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            cell_stresses[q_point] = C ⊡ ε
        end

        σ_avg = sum(cell_stresses) / getnquadpoints(cv)
        avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
        avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
        avg_cell_stresses[3][cellid(cell)] = σ_avg[1, 2]
    end

    return qp_stresses, avg_cell_stresses
end
# function calculate_stresses(grid, dh, cv, u, E, ν)
#     qp_stresses = [
#         [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
#         for _ in 1:getncells(grid)]
#     avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:3)...)
#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix(E[cell_index], ν)
#         cell_stresses = qp_stresses[cellid(cell)]
#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             cell_stresses[q_point] = C ⊡ ε
#         end
#         σ_avg = sum(cell_stresses) / getnquadpoints(cv)
#         avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
#         avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
#         avg_cell_stresses[3][cellid(cell)] = σ_avg[1, 2]
#     end
#     return qp_stresses, avg_cell_stresses
# end
# Function to calculate strains
"""
Function to calculate strains
```
calculate_strains(grid, dh, cv, u)
```
"""
function calculate_strains(grid, dh, cv, u)
    qp_strains = [
        [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
        for _ in 1:getncells(grid)]
    avg_cell_strains = tuple((zeros(getncells(grid)) for _ in 1:3)...)
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)
        cell_strains = qp_strains[cellid(cell)]
        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            cell_strains[q_point] = ε
        end
        ε_avg = sum(cell_strains) / getnquadpoints(cv)
        avg_cell_strains[1][cellid(cell)] = ε_avg[1, 1]
        avg_cell_strains[2][cellid(cell)] = ε_avg[2, 2]
        avg_cell_strains[3][cellid(cell)] = ε_avg[1, 2]
    end
    return qp_strains, avg_cell_strains
end
# Function to calculate element strain energy
"""
Function to calculate element strain energy
```
calculate_strain_energy(grid, dh, cv, u, E, ν)
```
"""
function calculate_strain_energy(grid, dh, cv, u, E, ν)
    element_strain_energies = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)
        
        # Check if E is scalar or vector
        local_E = isa(E, AbstractVector) ? E[cell_index] : E

        C = get_material_matrix(local_E, ν)
        cell_energy = 0.0
        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            σ = C ⊡ ε
            W = 0.5 * tr(σ ⊡ ε)
            dΩ = getdetJdV(cv, q_point)
            cell_energy += W * dΩ
        end
        element_strain_energies[cellid(cell)] = cell_energy
    end
    return element_strain_energies
end
# function calculate_strain_energy(grid, dh, cv, u, E, ν)
#     element_strain_energies = zeros(getncells(grid))
#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix(E[cell_index], ν)
#         cell_energy = 0.0
#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             σ = C ⊡ ε
#             W = 0.5 * tr(σ ⊡ ε)
#             dΩ = getdetJdV(cv, q_point)
#             cell_energy += W * dΩ
#         end
#         element_strain_energies[cellid(cell)] = cell_energy
#     end
#     return element_strain_energies
# end

"""
Function to calculate the derivative of the strain energy with respect to the Young's modulus
```
get_material_matrix_derivative_wrt_E(E, ν)
```
"""
function get_material_matrix_derivative_wrt_E(E, ν)
    dC_dE_voigt = [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, dC_dE_voigt)
end

# Function to calculate the cell volume
"""
Function to calculate the cell volume
```
calculate_cell_volume(cv)
```
"""
function calculate_cell_volume(cv)
    cell_volume = 0.0
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        cell_volume += dΩ
    end
    return cell_volume
end
# Function to calculate parameter H
"""
Function to calculate parameter H
```
calculate_H(grid, dh, cv, u, E, ν)
```
"""
function calculate_H(grid, dh, cv, u, E, ν)
    element_strain_energy_derivatives = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Handle E as either scalar or vector
        dC_dE = typeof(E) <: AbstractVector ? 
                get_material_matrix_derivative_wrt_E(E[cell_index], ν) : 
                get_material_matrix_derivative_wrt_E(E, ν)

        cell_derivative = 0.0
        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            dσ_dE = dC_dE ⊡ ε
            dW_dE = 0.5 * tr(dσ_dE ⊡ ε)
            dΩ = getdetJdV(cv, q_point)
            cell_derivative += dW_dE * dΩ
        end
        cell_volume = calculate_cell_volume(cv)
        element_strain_energy_derivatives[cellid(cell)] = cell_derivative / cell_volume
    end
    return element_strain_energy_derivatives
end
# function calculate_H(grid, dh, cv, u, E, ν)
#     element_strain_energy_derivatives = zeros(getncells(grid))
#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         dC_dE = get_material_matrix_derivative_wrt_E(E[cell_index], ν)
#         cell_derivative = 0.0
#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             dσ_dE = dC_dE ⊡ ε
#             dW_dE = 0.5 * tr(dσ_dE ⊡ ε)
#             dΩ = getdetJdV(cv, q_point)
#             cell_derivative += dW_dE * dΩ
#         end
#         cell_volume = calculate_cell_volume(cv)
#         element_strain_energy_derivatives[cellid(cell)] = cell_derivative / cell_volume
#     end
#     return element_strain_energy_derivatives
# end
######### average strain energy
"""
Function to calculate the average strain energy
```
calculate_average_strain_energy(grid, dh, cv, u, E, ν)
```
"""
function calculate_average_strain_energy(grid, dh, cv, u, E, ν)
    element_strain_energies = zeros(getncells(grid))
    element_volumes = zeros(getncells(grid))

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Handle E as either scalar or vector
        C = typeof(E) <: AbstractVector ? 
            get_material_matrix(E[cell_index], ν) : 
            get_material_matrix(E, ν)

        cell_energy = 0.0
        cell_volume = 0.0
        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            σ = C ⊡ ε
            W = 0.5 * tr(σ ⊡ ε)
            dΩ = getdetJdV(cv, q_point)
            cell_energy += W * dΩ
            cell_volume += dΩ
        end
        element_strain_energies[cellid(cell)] = cell_energy / cell_volume
        element_volumes[cellid(cell)] = cell_volume
    end
    return element_strain_energies
end
# function calculate_average_strain_energy(grid, dh, cv, u, E, ν)
#     element_strain_energies = zeros(getncells(grid))
#     element_volumes = zeros(getncells(grid)) 

#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix(E[cell_index], ν)
#         cell_energy = 0.0
#         cell_volume = 0.0
#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             σ = C ⊡ ε
#             W = 0.5 * tr(σ ⊡ ε)
#             dΩ = getdetJdV(cv, q_point)
#             cell_energy += W * dΩ
#             cell_volume += dΩ  # 
#         end
#         element_strain_energies[cellid(cell)] = cell_energy / cell_volume
#         element_volumes[cellid(cell)] = cell_volume
#     end
#     return element_strain_energies  
# end

struct LoadCondition
    load_type::String
    load_data::Union{Nothing,Vector{Float64},Float64}
end

# Struct for finite element solver
struct FEMSolver
    u::Vector{Float64}
    compliance::Float64
    σ::Any
    ε::Any
    U::Vector{Float64}
    H::Vector{Float64}
    ψ_avg::Vector{Float64}
end

# Main finite element solver
"""
Main finite element solver for two dimension
```
fem_solver(par::DynamicParams)
```
"""
function fem_solver(par::DynamicParams)

    # Unpack parameters
    grid = par.grid
    cell_values = par.cell_values
    facet_values = par.facet_values
    dh = par.dh
    ch = par.ch
    Neumann_bc = get(par, :Neumann_bc, nothing)  # Default to `nothing` if not provided
    E = par.E
    ν = par.ν
    loads = get(par, :loads, [])  # Default to an empty vector if not provided

    # Allocate and assemble global stiffness matrix
    K = allocate_matrix(dh)
    assemble_global!(K, dh, cell_values, E, ν)

    # Initialize external force vector
    f_ext = zeros(ndofs(dh))

    # Process each load condition and apply to Neumann_bc, if provided
    if loads !== nothing
        for load in loads
            if load.load_type == "traction_load"
                if Neumann_bc !== nothing
                    assemble_external_forces!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
                else
                    error("Neumann boundary conditions are required for traction loads.")
                end
            elseif load.load_type == "nodal_load"
                apply_nodal_force!(grid, Neumann_bc, load.load_data, f_ext)
            elseif load.load_type == "pressure_load"
                if Neumann_bc !== nothing
                    assemble_external_pressure!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
                else
                    error("Neumann boundary conditions are required for pressure loads.")
                end
            else
                error("Unknown load type: $(load.load_type)")
            end
        end
    end

    # Apply constraints and solve the system
    Ferrite.apply!(K, f_ext, ch)
    u = K \ f_ext  # Solve linear system
    c = 0.5 * dot(f_ext, u)

    # Calculate derived quantities
    _, σ = calculate_stresses(grid, dh, cell_values, u, E, ν)
    _, ε = calculate_strains(grid, dh, cell_values, u)
    U = calculate_strain_energy(grid, dh, cell_values, u, E, ν)
    H = calculate_H(grid, dh, cell_values, u, E, ν)
    ψ_avg = calculate_average_strain_energy(grid, dh, cell_values, u, E, ν)

    return FEMSolver(u, c, σ, ε, U, H, ψ_avg)
end

# function fem_solver(par::DynamicParams)

#     # grid, cell_values, facet_values, dh, ch, Neumann_bc, E, ν, loads::Vector{LoadCondition})
#     grid = par.grid
#     cell_values = par.cell_values
#     facet_values = par.facet_values
#     dh = par.dh
#     ch = par.ch
#     Neumann_bc = par.Neumann_bc
#     E = par.E
#     ν = par.ν
#     loads = par.loads
#     # Allocate and assemble global stiffness matrix
#     K = allocate_matrix(dh)
#     assemble_global!(K, dh, cell_values, E, ν)

#     # Initialize external force vector
#     f_ext = zeros(ndofs(dh))

#     # Process each load condition and apply to Neumann_bc
#     for load in loads
#         if load.load_type == "traction_load"
#             assemble_external_forces!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
#         elseif load.load_type == "nodal_load"
#             apply_nodal_force!(grid, Neumann_bc, load.load_data, f_ext)
#         elseif load.load_type == "pressure_load"
#             assemble_external_pressure!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
#         else
#             error("Unknown load type: $(load.load_type)")
#         end
#     end

#     # Apply constraints and solve the system
#     Ferrite.apply!(K, f_ext, ch)
#     u = K \ f_ext  # Solve linear system
#     c = 0.5 * dot(f_ext, u)

#     # Calculate derived quantities
#     _, σ = calculate_stresses(grid, dh, cell_values, u, E, ν)
#     _, ε = calculate_strains(grid, dh, cell_values, u)
#     U = calculate_strain_energy(grid, dh, cell_values, u, E, ν)
#     H = calculate_H(grid, dh, cell_values, u, E, ν)
#     ψ_avg = calculate_average_strain_energy(grid, dh, cell_values, u, E, ν)
#     return FEMSolver(u, c, σ, ε, U, H, ψ_avg)
# end


# # Example: Defining load conditions for Neumann_bc
# loads = [
#     LoadCondition("traction_load", [10.0, 0.0]),  # Traction load on Neumann_bc
#     LoadCondition("pressure_load", 5.0),         # Pressure load on Neumann_bc
#     LoadCondition("nodal_load", [0.0, -100.0])   # Nodal load on Neumann_bc
# ]

# # Call the solver with the load conditions
# result = fem_solver(grid, cell_values, facet_values, dh, ch, Neumann_bc; E=E, ν=ν, loads=loads)

###########################
# 3D functions
### finite elements code for 3d cases
"""
function to apply nodal forces to external force vector in 3D
```
apply_nodal_force_3d!(grid, node_set_name, load_vector, f_ext)
```
"""
function apply_nodal_force_3d!(grid, node_set_name, load_vector, f_ext)
    # Get the nodes in the specified nodeset
    node_set = getnodeset(grid, node_set_name)

    # Apply the load vector to each node in the nodeset
    for node_id in node_set
        f_ext[3*node_id-2] += load_vector[1]  
        f_ext[3*node_id-1] += load_vector[2]  
        f_ext[3*node_id] += load_vector[3]      
    end
end
##############################################
##############################################

"""
function to get the material matrix for 3D
```
get_material_matrix_3d(E, ν)
```
"""
function get_material_matrix_3d(E, ν)
    # Define the 3D material stiffness matrix in Voigt notation
    C_voigt = E / ((1 + ν) * (1 - 2 * ν)) * [
        1-ν ν ν 0 0 0;
        ν 1-ν ν 0 0 0;
        ν ν 1-ν 0 0 0;
        0 0 0 (1-2*ν)/2 0 0;
        0 0 0 0 (1-2*ν)/2 0;
        0 0 0 0 0 (1-2*ν)/2
    ]
    return fromvoigt(SymmetricTensor{4,3}, C_voigt)
end
##############################################
##############################################
"""
function to assemble the local stiffness matrix for 3D
```
assemble_cell_3d!(ke, cell_values, E, ν)
```
"""
function assemble_cell_3d!(ke, cell_values, E, ν)
    # Get the material matrix for 3D
    C = get_material_matrix_3d(E, ν)

    for qp in 1:getnquadpoints(cell_values)
        dΩ = getdetJdV(cell_values, qp)

        for i in 1:getnbasefunctions(cell_values)
            ∇Ni = shape_gradient(cell_values, qp, i)

            for j in 1:getnbasefunctions(cell_values)
                ∇δNj = shape_symmetric_gradient(cell_values, qp, j)
                ke[i, j] += (∇Ni ⊡ C ⊡ ∇δNj) * dΩ
            end
        end
    end
    return ke
end
##############################################
##############################################
"""
function to assemble the global stiffness matrix for 3D
```
assemble_global_3d!(K, dh, cell_values, E, ν)
```
"""
function assemble_global_3d!(K, dh, cell_values, E, ν)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)

        # Handle E as either scalar or vector
        material_E = typeof(E) <: AbstractVector ? E[cell_index] : E
        assemble_cell_3d!(ke, cell_values, material_E, ν)

        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end
# function assemble_global_3d!(K, dh, cell_values, E, ν)
#     n_basefuncs = getnbasefunctions(cell_values)
#     ke = zeros(n_basefuncs, n_basefuncs)
#     assembler = start_assemble(K)

#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cell_values, cell)
#         fill!(ke, 0.0)
#         assemble_cell_3d!(ke, cell_values, E[cell_index], ν)
#         assemble!(assembler, celldofs(cell), ke)
#     end
#     return K
# end
########## traction surface forces
"""
function to assemble external forces from surface tractions in 3D
```
assemble_external_forces_3d!(f_ext, dh, facetset, facet_values, traction)
```
"""
function assemble_external_forces_3d!(f_ext, dh, facetset, facet_values, traction)
    fe_ext = zeros(getnbasefunctions(facet_values))
    for facet in FacetIterator(dh, facetset)
        reinit!(facet_values, facet)
        fill!(fe_ext, 0.0)
        for qp in 1:getnquadpoints(facet_values)
            dΓ = getdetJdV(facet_values, qp)
            for i in 1:getnbasefunctions(facet_values)
                fe_ext[i] += traction ⋅ shape_value(facet_values, qp, i) * dΓ
            end
        end
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
##### pressure surface forces
"""
function to assemble external forces from pressure in 3D
```
assemble_external_pressure_3d!(f_ext, dh, facetset, facet_values, pressure)
```
"""
function assemble_external_pressure_3d!(f_ext, dh, facetset, facet_values, pressure)
    fe_ext = zeros(getnbasefunctions(facet_values))
    for facet in FacetIterator(dh, facetset)
        reinit!(facet_values, facet)
        fill!(fe_ext, 0.0)
        for qp in 1:getnquadpoints(facet_values)
            p = pressure * getnormal(facet_values, qp)
            dΓ = getdetJdV(facet_values, qp)
            for i in 1:getnbasefunctions(facet_values)
                fe_ext[i] += p ⋅ shape_value(facet_values, qp, i) * dΓ
            end
        end
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
##### 
"""
function to calculate stresses in 3D
```
calculate_stresses_3d(grid, dh, cv, u, E, ν)
```
"""
function calculate_stresses_3d(grid, dh, cv, u, E, ν)
    # Initialize containers for stresses
    qp_stresses = [
        [zero(SymmetricTensor{2,3}) for _ in 1:getnquadpoints(cv)]
        for _ in 1:getncells(grid)]
    avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:6)...)  # 3D Voigt format

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Handle E as either scalar or vector
        C = typeof(E) <: AbstractVector ? 
            get_material_matrix_3d(E[cell_index], ν) : 
            get_material_matrix_3d(E, ν)

        cell_stresses = qp_stresses[cellid(cell)]

        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            cell_stresses[q_point] = C ⊡ ε
        end

        # Average stresses for the cell
        σ_avg = sum(cell_stresses) / getnquadpoints(cv)
        avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
        avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
        avg_cell_stresses[3][cellid(cell)] = σ_avg[3, 3]
        avg_cell_stresses[4][cellid(cell)] = σ_avg[1, 2]
        avg_cell_stresses[5][cellid(cell)] = σ_avg[2, 3]
        avg_cell_stresses[6][cellid(cell)] = σ_avg[1, 3]
    end
    return qp_stresses, avg_cell_stresses
end
# function calculate_stresses_3d(grid, dh, cv, u, E, ν)
#     # Initialize containers for stresses
#     qp_stresses = [
#         [zero(SymmetricTensor{2,3}) for _ in 1:getnquadpoints(cv)]
#         for _ in 1:getncells(grid)]
#     avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:6)...)  # 3D Voigt format

#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix_3d(E[cell_index], ν)
#         cell_stresses = qp_stresses[cellid(cell)]

#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             cell_stresses[q_point] = C ⊡ ε
#         end

#         # Average stresses for the cell
#         σ_avg = sum(cell_stresses) / getnquadpoints(cv)
#         avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
#         avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
#         avg_cell_stresses[3][cellid(cell)] = σ_avg[3, 3]
#         avg_cell_stresses[4][cellid(cell)] = σ_avg[1, 2]
#         avg_cell_stresses[5][cellid(cell)] = σ_avg[2, 3]
#         avg_cell_stresses[6][cellid(cell)] = σ_avg[1, 3]
#     end
#     return qp_stresses, avg_cell_stresses
# end
########## strains
##############################################################################
##############################################################################
##############################################################################
"""
function to calculate strains in 3D
```
calculate_strains_3d(grid, dh, cv, u)
```
"""
function calculate_strains_3d(grid, dh, cv, u)
    # Initialize containers for strains
    qp_strains = [
        [zero(SymmetricTensor{2,3}) for _ in 1:getnquadpoints(cv)]
        for _ in 1:getncells(grid)]
    avg_cell_strains = tuple((zeros(getncells(grid)) for _ in 1:6)...)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)
        cell_strains = qp_strains[cellid(cell)]

        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            cell_strains[q_point] = ε
        end

        # Average strains for the cell
        ε_avg = sum(cell_strains) / getnquadpoints(cv)
        avg_cell_strains[1][cellid(cell)] = ε_avg[1, 1]
        avg_cell_strains[2][cellid(cell)] = ε_avg[2, 2]
        avg_cell_strains[3][cellid(cell)] = ε_avg[3, 3]
        avg_cell_strains[4][cellid(cell)] = ε_avg[1, 2]
        avg_cell_strains[5][cellid(cell)] = ε_avg[2, 3]
        avg_cell_strains[6][cellid(cell)] = ε_avg[1, 3]
    end
    return qp_strains, avg_cell_strains
end

##############################################################################
##############################################################################
##############################################################################
"""
function to calculate element strain energy in 3D
```
calculate_strain_energy_3d(grid, dh, cv, u, E, ν)
```
"""
function calculate_strain_energy_3d(grid, dh, cv, u, E, ν)
    # Initialize the strain energy array
    element_strain_energies = zeros(getncells(grid))

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Handle E as either scalar or vector
        C = typeof(E) <: AbstractVector ? 
            get_material_matrix_3d(E[cell_index], ν) : 
            get_material_matrix_3d(E, ν)

        cell_energy = 0.0

        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            σ = C ⊡ ε
            W = 0.5 * tr(σ ⊡ ε)
            dΩ = getdetJdV(cv, q_point)
            cell_energy += W * dΩ
        end

        element_strain_energies[cellid(cell)] = cell_energy
    end
    return element_strain_energies
end
# function calculate_strain_energy_3d(grid, dh, cv, u, E, ν)
#     # Initialize the strain energy array
#     element_strain_energies = zeros(getncells(grid))

#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix_3d(E[cell_index], ν)
#         cell_energy = 0.0

#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             σ = C ⊡ ε
#             W = 0.5 * tr(σ ⊡ ε)
#             dΩ = getdetJdV(cv, q_point)
#             cell_energy += W * dΩ
#         end

#         element_strain_energies[cellid(cell)] = cell_energy
#     end
#     return element_strain_energies
# end

"""
function to get the derivative of the material matrix with respect to the Young's modulus in 3D
```
get_material_matrix_derivative_wrt_E_3d(E::Vector{Float64,1}, ν)
```
"""
function get_material_matrix_derivative_wrt_E_3d(E, ν)
    dC_dE_voigt = [
        1.0 ν ν 0.0 0.0 0.0;
        ν 1.0 ν 0.0 0.0 0.0;
        ν ν 1.0 0.0 0.0 0.0;
        0.0 0.0 0.0 (1-ν)/2.0 0.0 0.0;
        0.0 0.0 0.0 0.0 (1-ν)/2.0 0.0;
        0.0 0.0 0.0 0.0 0.0 (1-ν)/2.0
    ] / (1 - 2 * ν) / (1 + ν)
    return fromvoigt(SymmetricTensor{4,3}, dC_dE_voigt)
end
function calculate_cell_volume_3d(cv)
    cell_volume = 0.0
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        cell_volume += dΩ
    end
    return cell_volume
end

"""
function to calculate parameter H in 3D
```
calculate_H_3d(grid, dh, cv, u, E, ν)
```
"""
function calculate_H_3d(grid, dh, cv, u, E, ν)
    # Initialize the derivatives array
    element_strain_energy_derivatives = zeros(getncells(grid))

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Determine the elasticity modulus for the current cell
        E_cell = isa(E, Number) ? E : E[cell_index]

        # Get the material matrix and its derivative for the current cell
        C = get_material_matrix_3d(E_cell, ν)
        dC_dE = get_material_matrix_derivative_wrt_E_3d(E_cell, ν)

        cell_derivative = 0.0

        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            dσ_dE = dC_dE ⊡ ε  # ⊗ represents the double contraction operation
            dW_dE = 0.5 * tr(dσ_dE ⊡ ε)
            dΩ = getdetJdV(cv, q_point)
            cell_derivative += dW_dE * dΩ
        end

        cell_volume = calculate_cell_volume_3d(cv)
        element_strain_energy_derivatives[cellid(cell)] = cell_derivative / cell_volume
    end

    return element_strain_energy_derivatives
end
# function calculate_H_3d(grid, dh, cv, u, E, ν)
#     # Initialize the derivatives array
#     element_strain_energy_derivatives = zeros(getncells(grid))

#     for (cell_index, cell) in enumerate(CellIterator(dh))
#         reinit!(cv, cell)
#         C = get_material_matrix_3d(E[cell_index], ν)
#         dC_dE = get_material_matrix_derivative_wrt_E_3d(E[cell_index], ν)
#         cell_derivative = 0.0

#         for q_point in 1:getnquadpoints(cv)
#             ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
#             dσ_dE = dC_dE ⊡ ε
#             dW_dE = 0.5 * tr(dσ_dE ⊡ ε)
#             dΩ = getdetJdV(cv, q_point)
#             cell_derivative += dW_dE * dΩ
#         end

#         cell_volume = calculate_cell_volume_3d(cv)
#         element_strain_energy_derivatives[cellid(cell)] = cell_derivative / cell_volume
#     end
#     return element_strain_energy_derivatives
# end
# Struct for load conditions in 3D
struct LoadCondition_3d
    load_type::String
    load_data::Union{Nothing,Vector{Float64},Float64}
end

# Struct for finite element solver results in 3D
struct FEMSolver_3d
    u::Vector{Float64}        # Displacement vector
    compliance::Float64       # Compliance value
    σ::Any                    # Stresses (can be a tensor or another representation)
    ε::Any                    # Strains (can be a tensor or another representation)
    U::Vector{Float64}        # Strain energy vector
    H::Vector{Float64}        # Other derived quantity (e.g., sensitivity or residual)
end

# Main finite element solver
"""
Main finite element solver for 3D
```
fem_solver_3d(par::DynamicParams)
```
"""
function fem_solver_3d(par::DynamicParams)
    # Unpack parameters
    grid = par.grid
    cell_values = par.cell_values
    facet_values = par.facet_values
    dh = par.dh
    ch = par.ch
    Neumann_bc = get(par, :Neumann_bc, nothing)  # Default to `nothing` if not provided
    E = par.E
    ν = par.ν
    loads = get(par, :loads, [])  # Default to an empty vector if not provided

    # Allocate and assemble global stiffness matrix
    K = allocate_matrix(dh)  # Allocate stiffness matrix
    assemble_global_3d!(K, dh, cell_values, E, ν)

    # Initialize external force vector
    f_ext = zeros(ndofs(dh))

    # Process each load condition and apply to Neumann boundary conditions, if provided
    if loads !== nothing
        for load in loads
            if load.load_type == "traction_load"
                if Neumann_bc !== nothing
                    assemble_external_forces_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
                else
                    error("Neumann boundary conditions are required for traction loads.")
                end
            elseif load.load_type == "nodal_load"
                apply_nodal_force_3d!(grid, Neumann_bc, load.load_data, f_ext)
            elseif load.load_type == "pressure_load"
                if Neumann_bc !== nothing
                    assemble_external_pressure_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
                else
                    error("Neumann boundary conditions are required for pressure loads.")
                end
            else
                error("Unknown load type: $(load.load_type)")
            end
        end
    end

    # Apply constraints and solve the linear system
    apply!(K, f_ext, ch)  # Apply constraints
    u = K \ f_ext         # Solve the system Ku = f_ext
    c = 0.5 * dot(f_ext, u)  # Compute compliance

    # Calculate derived quantities
    _, σ = calculate_stresses_3d(grid, dh, cell_values, u, E, ν)  # Stresses
    _, ε = calculate_strains_3d(grid, dh, cell_values, u)         # Strains
    U = calculate_strain_energy_3d(grid, dh, cell_values, u, E, ν)  # Strain energy
    H = calculate_H_3d(grid, dh, cell_values, u, E, ν)            # Derived quantity H

    # Return solver results
    return FEMSolver_3d(u, c, σ, ε, U, H)
end

# function fem_solver_3d(par::DynamicParams)
#     #grid, cell_values, facet_values, dh, ch, Neumann_bc, E, ν, loads::Vector{LoadCondition_3d}
#     grid = par.grid
#     cell_values = par.cell_values
#     facet_values = par.facet_values
#     dh = par.dh
#     ch = par.ch
#     Neumann_bc = par.Neumann_bc
#     E = par.E
#     ν = par.ν
#     loads = par.loads
#     # Allocate and assemble global stiffness matrix
#     K = allocate_matrix(dh)  # Allocate stiffness matrix
#     assemble_global_3d!(K, dh, cell_values, E, ν)

#     # Initialize external force vector
#     f_ext = zeros(ndofs(dh))

#     # Process each load condition and apply to Neumann boundary conditions
#     for load in loads
#         if load.load_type == "traction_load"
#             assemble_external_forces_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
#         elseif load.load_type == "nodal_load"
#             apply_nodal_force_3d!(grid, Neumann_bc, load.load_data, f_ext)
#         elseif load.load_type == "pressure_load"
#             assemble_external_pressure_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
#         else
#             error("Unknown load type: $(load.load_type)")
#         end
#     end

#     # Apply constraints and solve the linear system
#     apply!(K, f_ext, ch)  # Apply constraints
#     u = K \ f_ext         # Solve the system Ku = f_ext
#     c = 0.5 * dot(f_ext, u)  # Compute compliance

#     # # Calculate derived quantities
#      _, σ = calculate_stresses_3d(grid, dh, cell_values, u, E, ν)  # Stresses
#     _, ε = calculate_strains_3d(grid, dh, cell_values, u)         # Strains
#     U = calculate_strain_energy_3d(grid, dh, cell_values, u, E, ν)  # Strain energy
#     H = calculate_H_3d(grid, dh, cell_values, u, E, ν)            # Derived quantity H

#     # Return solver results
#     return FEMSolver_3d(u, c, σ, ε, U, H)
# end

###################### functions for topology optimization based on UPM approach


"""
function for update Young's modulus based on UPM approach

Example:
````
k = 1
E = rand(Float64, 5)
H = rand(Float64, 5)
Emax = 1.0
Emin = 1e-4
# test code
Enew = update_upm!(k, E, H,Emax,Emin)

````
"""
function update_upm!(k::Int64, E::Array{Float64,1}, H::Array{Float64,1}, Emax::Float64, Emin::Float64)
    Enew = copy(E)
    # H_mean = calculate_mean(H)
    # H_std = standard_deviation(H)
    H_mean = mean(H)
    H_std = std(H)
    for i in eachindex(E)
        α = (H[i] - H_mean) / (k * H_std)
        Enew[i] = E[i] * (1 + α)
        Enew = [clamp(x, Emin, Emax) for x in Enew]
    end
    return Enew
end

"""
function to convert Young's modulus to density
    E = E0(rho/rho0)^gamma

Example:
````
Enew = rand(Float64, 5)
E0 = 1.0
ρ0 = 1.0
γ = 1
ρ = transfer_to_density!(Enew, E0, ρ0, γ)


````
"""
function transfer_to_density!(Enew::Array{Float64,1}, E0::Float64, ρ0::Float64, γ::Int64)
    ρmin, ρmax = 0.0, 1.0

    ρ = ρ0 .* ((Enew / E0) .^ (1 / γ))

    ρ = [clamp(x, ρmin, ρmax) for x in ρ]

    return ρ
end

"""
function to transfer density to young modulus
example:
ρnew = rand(Float64, 4)
E0 = 1.0
ρ0 = 1.0
γ = 1
Emin = 1e-4
Emax = 1.0

transfer_to_young!(ρnew , E0, ρ0, γ, Emin , Emax)

"""
function transfer_to_young!(ρnew::Array{Float64,1}, E0::Float64,
    ρ0::Float64, γ::Int64, Emin::Float64, Emax::Float64)

    Enew = E0 * (ρnew / ρ0) .^ γ
    Enew = [clamp(x, Emin, Emax) for x in Enew]

    return Enew
end

"""
function to apply volume fraction

example:
```
ρnew = rand(Float64, 8)
nx, ny , nz = 2 , 2 , 2
vf = 0.5
η = π/4

ρ =  filter_density_to_vf!(ρnew, vf, nx, ny, nz, η)
```
"""
function filter_density_to_vf!(density, vf, nx, ny, nz, eta)
    rhomin, rhomax = 0.01, 1.
    function transform(rholoc, rhotr, eta, rhomin, rhomax)
        if rholoc < rhotr
            rhotrans = rhomin  
        elseif rholoc > rhotr + 1.0/tan(eta)
            rhotrans = rhomax
        else
            rhotrans = tan(eta) * (rholoc - rhotr)
        end
        return rhotrans
    end
    rhomaxbound = -1.0/tan(eta)  # minimum that gives a vf of 0
    rhominbound = 1.0  # maximum that gives a vf of 1
    error = 10.0  # just put a high number
    rhotr = 0.0  # Initialize rhotr before the loop
    while error > 0.001
        rhotr = (rhominbound + rhomaxbound) / 2  # this is the initial point
        sumdmin = 0.0
        for i in eachindex(density)
            sumdmin += transform(density[i], rhominbound, eta, rhomin, rhomax) / (nx * ny * nz)
        end
        sumdmax = 0.0
        for i in eachindex(density)
            sumdmax += transform(density[i], rhomaxbound, eta, rhomin, rhomax) / (nx * ny * nz)
        end
        sumdmid = 0.0
        for i in eachindex(density)
            sumdmid += transform(density[i], rhotr, eta, rhomin, rhomax) / (nx * ny * nz)
        end
        if (sumdmin - vf) / (sumdmid - vf) > 0
            rhominbound = rhotr
        elseif (sumdmax - vf) / (sumdmid - vf) > 0
            rhomaxbound = rhotr
        else
            println("problem out of bounds", sumdmax, sumdmin, sumdmid, vf)
        end
        error = abs(vf - sumdmid)
    end
    for i in eachindex(density)
        densloc = transform(density[i], rhotr, eta, rhomin, rhomax)
        density[i] = densloc
    end
    return density
end

"""
function to perform topology optimization using UPM approach (2D case)
```
top_upm!(par::DynamicParams, name_of_file::String, directory::String)
```    
"""
function top_upm!(par::DynamicParams, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    E = par.E
    nx = par.nx ; ny = par.ny ; nz = par.nz
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    k = par.k ; γ = par.γ ; volfrac = par.vf; η = par.η; ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver(par)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

    println("Iter $loop: C = $compliance, ΔR = $A")

    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        Enew = update_upm!(k, E, H, Emax, Emin)
        ρ = transfer_to_density!(Enew, E0, ρ0, γ)
        ρnew = filter_density_to_vf!(ρ, volfrac, nx, ny, nz, η)
        Enew_frac = transfer_to_young!(ρnew, E0, ρ0, γ, Emin, Emax)

        # Update E in par so that fem_solver uses the updated material distribution
        par.E = Enew_frac
        E = Enew_frac

        fem = fem_solver(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        u = fem.u
        σ = fem.σ
        ε = fem.ε

        # Create iteration-specific filename with zero-padded loop index
        loop_str = @sprintf("%03d", loop)
        filename_it = string(name_of_file, "_", loop_str)
        full_path = joinpath(directory, filename_it)

        VTKGridFile(full_path, dh) do vtk
            write_solution(vtk, dh, u)
            for (j, key) in enumerate(("11", "22", "12"))
                write_cell_data(vtk, σ[j], "stress_" * key)
            end
            for (j, key) in enumerate(("11", "22", "12"))
                write_cell_data(vtk, ε[j], "strain_" * key)
            end
            write_cell_data(vtk, E, "Young's modulus")
            Ferrite.write_cellset(vtk, grid)
        end

        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")

    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver(par)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        for (j, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, par.E, "Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end

end
    
"""
function to perform topology optimization using UPM approach (3D case)
```
top_upm_3d!(par::DynamicParams, name_of_file::String, directory::String)
```
"""
function top_upm_3d!(par::DynamicParams, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    E = par.E
    nx = par.nx ; ny = par.ny ; nz = par.nz
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    k = par.k ; γ = par.γ ; volfrac = par.vf ; η = par.η ; ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    # Initial FEM solve
    fem = fem_solver_3d(par)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)

    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")

    # Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        # FEM solve with current parameters
        fem = fem_solver_3d(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        # Material update routines
        Enew = update_upm!(k, E, H, Emax, Emin)
        ρ = transfer_to_density!(Enew, E0, ρ0, γ)
        ρnew = filter_density_to_vf!(ρ, volfrac, nx, ny, nz, η)
        Enew_frac = transfer_to_young!(ρnew, E0, ρ0, γ, Emin, Emax)

        # Update par to reflect the new E distribution
        par.E = Enew_frac
        E = Enew_frac

        # FEM solve after material update
        fem = fem_solver_3d(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        u = fem.u
        σ = fem.σ
        ε = fem.ε

        # Create iteration-specific filename
        loop_str = @sprintf("%03d", loop)
        filename_it = string(name_of_file, "_", loop_str)
        full_path = joinpath(directory, filename_it)

        VTKGridFile(full_path, dh) do vtk
            write_solution(vtk, dh, u)
            # Write 3D stress components
            for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
                write_cell_data(vtk, σ[j], "stress_" * key)
            end
            # Write 3D strain components
            for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
                write_cell_data(vtk, ε[j], "strain_" * key)
            end
            write_cell_data(vtk, E, "Young's modulus")
            Ferrite.write_cellset(vtk, grid)
        end

        # Update strain energy vector and A
        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")
    end

    # Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final solve and write results
    fem = fem_solver_3d(par)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        # Write 3D stress components
        for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        # Write 3D strain components
        for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end
end