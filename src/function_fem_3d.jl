############################################

# 3D functions
### finite elements code for 3d cases
"""
function to apply nodal forces to external force vector in 3D
```
apply_nodal_force_3d!(grid, node_set_name, load_vector, f_ext)
```
"""
function apply_nodal_force_3d!(grid, node_ids, load_vector, f, dh)
    # Get the coordinates of the nodes
    node_coords = [Ferrite.get_node_coordinate(grid, id) for id in node_ids]

    # Detect the primary edge direction by checking which coordinate is approximately constant
    tolerance = 1e-8
    is_x_constant = all(abs(node_coords[1][1] - coord[1]) < tolerance for coord in node_coords)
    is_y_constant = all(abs(node_coords[1][2] - coord[2]) < tolerance for coord in node_coords)
    is_z_constant = all(abs(node_coords[1][3] - coord[3]) < tolerance for coord in node_coords)

    # Determine the primary and secondary directions
    if is_x_constant
        primary_direction = 2  # y-direction
        secondary_direction = 3  # z-direction
    elseif is_y_constant
        primary_direction = 1  # x-direction
        secondary_direction = 3  # z-direction
    elseif is_z_constant
        primary_direction = 1  # x-direction
        secondary_direction = 2  # y-direction
    else
        error("The edge is not aligned with a primary axis!")
    end

    # Convert node IDs (OrderedSet) to a Vector for indexing
    node_id_vec = collect(node_ids)

    # Handle the case for a single node
    if length(node_ids) == 1
        # Get the single node ID
        single_node = node_ids[1]

        # Get the vertex index for this node
        vertex = nodeid_to_vertexindex(grid, single_node)

        # Get the DOFs associated with this vertex
        dofs = vertexdofs(dh, vertex)

        # Apply the entire load vector directly to this node
        for j in 1:3  # Iterate over x, y, z directions
            f[dofs[j]] += load_vector[j]
        end
    else
        # Sort nodes based on the primary and secondary directions for consistent ordering
        sorted_indices = sortperm([(coord[primary_direction], coord[secondary_direction]) for coord in node_coords])
        sorted_node_ids = node_id_vec[sorted_indices]

        # Get the sorted coordinates for the nodes
        sorted_coords = [node_coords[idx] for idx in sorted_indices]

        # Compute segment lengths in full 3D space
        segment_lengths = [
            norm([
                sorted_coords[i+1][j] - sorted_coords[i][j] for j in 1:3
            ]) for i in 1:(length(sorted_coords)-1)
        ]

        # Map forces to the DOFs of nodes
        total_length = sum(segment_lengths)
        for i in eachindex(segment_lengths)
            # Get the node IDs for the current segment
            node1 = sorted_node_ids[i]
            node2 = sorted_node_ids[i+1]

            # Calculate the segment's contribution to the total force
            segment_force = load_vector .* (segment_lengths[i] / total_length)  # Proportional to segment length

            # Get the vertex indices for these nodes
            vertex1 = nodeid_to_vertexindex(grid, node1)
            vertex2 = nodeid_to_vertexindex(grid, node2)

            # Get the DOFs associated with these vertices
            dofs1 = vertexdofs(dh, vertex1)
            dofs2 = vertexdofs(dh, vertex2)

            # Distribute forces equally to the DOFs of node1 and node2
            for j in 1:3  # Iterate over x, y, z directions
                f[dofs1[j]] += segment_force[j] / 2  # Half the force to node1
                f[dofs2[j]] += segment_force[j] / 2  # Half the force to node2
            end
        end
    end

    return f
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
function compute_nodal_data_3D(grid, element_data)
    nnodes = Ferrite.getnnodes(grid)  # Total number of nodes
    
    nodal_data = zeros(Float64, nnodes)  # Initialize nodal data array
    count_elements = zeros(Int, nnodes)  # To count elements connected to each node

    cells = getcells(grid)  # Get all cells
    for (element_id, cell) in enumerate(cells)
        node_ids = cell.nodes  # Access node indices directly from the cell
        for node_id in node_ids
            nodal_data[node_id] += element_data[element_id]
            count_elements[node_id] += 1
        end
    end

    # Average the sum of element data by the number of connected elements
    for node_id in 1:nnodes
        if count_elements[node_id] > 0
            nodal_data[node_id] /= count_elements[node_id]
        end
    end

    return nodal_data
end

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
    E_node::Vector{Float64}
end

# Main finite element solver
"""
Main finite element solver for 3D
```
fem_solver_3d(par::DynamicParams)
```
"""
function fem_solver_3d(par::DynamicParams)
    # Extract common parameters from the input
    grid = par.grid
    cell_values = par.cell_values
    facet_values = par.facet_values
    dh = par.dh
    ch = par.ch
    E = par.E
    ν = par.ν

    # Allocate and assemble the global stiffness matrix
    K = allocate_matrix(dh)
    assemble_global_3d!(K, dh, cell_values, E, ν)

    # Initialize external force vector
    f_ext = zeros(ndofs(dh))

    # Check if Neumann_bc and loads are not empty
    if !isempty(par.Neumann_bc) && !isempty(par.loads)
        Neumann_bc = par.Neumann_bc
        loads = par.loads

        for load in loads
            if load.load_type == "traction_load"
                assemble_external_forces_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
            elseif load.load_type == "nodal_load"
                apply_nodal_force_3d!(grid, Neumann_bc, load.load_data, f_ext, dh)
            elseif load.load_type == "pressure_load"
                assemble_external_pressure_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
            else
                error("Unknown load type: $(load.load_type)")
            end
        end
    end

    # Handle the case when Neumann_bc or loads are empty
    if isempty(par.Neumann_bc) || isempty(par.loads)
        println("Warning: Either Neumann boundary conditions or loads are empty. Proceeding without external forces.")
    end

    # Apply constraints and solve the system
    apply!(K, f_ext, ch)  # Apply constraints
    u = K \ f_ext         # Solve the system Ku = f_ext
    c = 0.5 * dot(f_ext, u)  # Compute compliance

    # Calculate derived quantities
    _, σ = calculate_stresses_3d(grid, dh, cell_values, u, E, ν)  # Stresses
    _, ε = calculate_strains_3d(grid, dh, cell_values, u)         # Strains
    U = calculate_strain_energy_3d(grid, dh, cell_values, u, E, ν)  # Strain energy
    H = calculate_H_3d(grid, dh, cell_values, u, E, ν)            # Derived quantity H
    
    E_node = compute_nodal_data_3D(grid, E)

    # Return solver results
    return FEMSolver_3d(u, c, σ, ε, U, H,E_node)
end