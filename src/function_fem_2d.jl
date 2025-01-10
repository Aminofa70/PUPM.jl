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
function vertexdofs(dh::DofHandler, vertexid::VertexIndex)
    cellid, lvidx = vertexid
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cellid]]
    local_vertex_dofs = Int[]

    for ifield in 1:length(sdh.field_names)
        offset = Ferrite.field_offset(sdh, ifield)
        field_dim = Ferrite.n_components(sdh, ifield)
        field_ip = isa(sdh.field_interpolations[ifield], Ferrite.VectorizedInterpolation) ?
                   sdh.field_interpolations[ifield].ip :
                   sdh.field_interpolations[ifield]

        vert = Ferrite.vertexdof_indices(field_ip)[lvidx]

        for vdof in vert, d in 1:field_dim
            push!(local_vertex_dofs, (vdof - 1) * field_dim + d + offset)
        end
    end

    dofs = zeros(Int, ndofs_per_cell(dh, cellid))
    celldofs!(dofs, dh, cellid)

    return dofs[local_vertex_dofs]
end
############################################
function nodeid_to_vertexindex(grid::Grid, nodeid::Int)
    for (cellid, cell) in enumerate(grid.cells)
        for (i, nodeid2) in enumerate(cell.nodes)
            if nodeid == nodeid2
                return VertexIndex(cellid, i)
            end
        end
    end
    error("Node $(nodeid) does not belong to any cell")
end
############################################
function apply_nodal_force!(grid, nodeid, load_vector, f, dh)
    # Get the coordinates of the nodes
    coords = [Ferrite.get_node_coordinate(grid, id) for id in nodeid]

    # Determine if the edge is vertical (x constant) or horizontal (y constant)
    is_vertical = all(abs(coords[1][1] - coord[1]) < 1e-8 for coord in coords)

    # Set primary direction based on edge orientation
    primary_direction = is_vertical ? 2 : 1  # 2 = y-direction, 1 = x-direction

    # Convert OrderedSet to Vector for indexing
    nodeid_vec = collect(nodeid)

    # Check if there is only one node
    if length(nodeid) == 1
        # Handle the single-node case
        single_node = nodeid[1]

        # Get the vertex index and DOFs for this node
        vertex = nodeid_to_vertexindex(grid, single_node)
        dofs = vertexdofs(dh, vertex)

        # Apply the entire load vector directly to this node's DOFs
        f[dofs[1:2]] .+= load_vector
    else
        # Handle the multi-node case
        # Sort nodes based on the primary direction
        sorted_indices = sortperm([coord[primary_direction] for coord in coords])
        sorted_nodeid = nodeid_vec[sorted_indices]

        # Compute segment lengths (dy or dx) along the primary direction
        sorted_coords = [coords[idx] for idx in sorted_indices]
        segment_lengths = diff([coord[primary_direction] for coord in sorted_coords])

        # Calculate the total length of the edge for proportional force distribution
        total_length = sum(segment_lengths)

        # Distribute forces among nodes proportionally to segment lengths
        for i in eachindex(segment_lengths)
            # Get the node IDs for the current segment
            node1 = sorted_nodeid[i]
            node2 = sorted_nodeid[i + 1]

            # Calculate the segment's proportional contribution to the load vector
            segment_force = load_vector * segment_lengths[i] / total_length

            # Get the vertex indices and DOFs for the nodes
            vertex1 = nodeid_to_vertexindex(grid, node1)
            vertex2 = nodeid_to_vertexindex(grid, node2)
            dofs1 = vertexdofs(dh, vertex1)
            dofs2 = vertexdofs(dh, vertex2)

            # Distribute half of the segment's force to each node
            f[dofs1[1:2]] .+= segment_force / 2
            f[dofs2[1:2]] .+= segment_force / 2
        end
    end

    return f
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
function compute_nodal_data(grid, element_data)
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
    E_node::Vector{Float64}
end

# Main finite element solver
"""
Main finite element solver for two dimension
```
fem_solver(par::DynamicParams)
```
"""
function fem_solver(par::DynamicParams)
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
    assemble_global!(K, dh, cell_values, E, ν)

    # Initialize external force vector
    f_ext = zeros(ndofs(dh))

    # Check if Neumann_bc and loads are not empty
    if !isempty(par.Neumann_bc) && !isempty(par.loads)
        Neumann_bc = par.Neumann_bc
        loads = par.loads

        for load in loads
            if load.load_type == "traction_load"
                assemble_external_forces!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
            elseif load.load_type == "nodal_load"
                apply_nodal_force!(grid, Neumann_bc, load.load_data, f_ext, dh)
            elseif load.load_type == "pressure_load"
                assemble_external_pressure!(f_ext, dh, Neumann_bc, facet_values, load.load_data)
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
    Ferrite.apply!(K, f_ext, ch)
    u = K \ f_ext  # Solve the linear system
    c = 0.5 * dot(f_ext, u)

    # Calculate derived quantities
    _, σ = calculate_stresses(grid, dh, cell_values, u, E, ν)
    _, ε = calculate_strains(grid, dh, cell_values, u)
    U = calculate_strain_energy(grid, dh, cell_values, u, E, ν)
    H = calculate_H(grid, dh, cell_values, u, E, ν)
    ψ_avg = calculate_average_strain_energy(grid, dh, cell_values, u, E, ν)
    E_node = compute_nodal_data(grid, E)


    # Return the result
    return FEMSolver(u, c, σ, ε, U, H, ψ_avg, E_node)
end
##########################################################
# # Example: Defining load conditions for Neumann_bc
# loads = [
#     LoadCondition("traction_load", [10.0, 0.0]),  # Traction load on Neumann_bc
#     LoadCondition("pressure_load", 5.0),         # Pressure load on Neumann_bc
#     LoadCondition("nodal_load", [0.0, -100.0])   # Nodal load on Neumann_bc
# ]

# # Call the solver with the load conditions
# result = fem_solver(grid, cell_values, facet_values, dh, ch, Neumann_bc; E=E, ν=ν, loads=loads)
##########################################################
###########################