using Revise
using Ferrite
using PUPM

## cantilever for UPM
par = DynamicParams()

# Function to create a 3D grid
function create_grid(L, w, h, nx, ny, nz, r = 0.1)
    # Define the domain boundaries using Vec{3}
    left = Vec(0.0, 0.0, 0.0)
    right = Vec(L, w, h)
    # Generate the 3D hexahedral grid
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    # Define a small load region in the center of the y-z plane at x = L
    addfacetset!(grid, "load", x -> x[1] ≈ L && abs(x[2] - w/2) <= r && abs(x[3] - h/2) <= r)
    return grid
end

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron, order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    dbcs = ConstraintHandler(dh)
    # Clamped on the left side
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], dofs)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end;

# Main script
L, w, h = 6.0, 1.0, 1.0
nx, ny, nz = 35, 10, 10

grid = create_grid(L, w, h, nx, ny, nz)

par.tnele = length(grid.cells)  # Total number of elements

par.grid = grid
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh, grid)

par.cell_values, par.facet_values = create_values()

#par.loads = [LoadCondition_3d("nodal_load", [0.0, 0.0, 1.0])]
par.loads = [LoadCondition_3d("traction_load", [0.0, 0.0, -1.0])]


# Material properties
par.E0 = 1.0
E = fill(par.E0,Ferrite.getncells(grid))
par.ν = 0.3

# Optimization parameters
par.Emin = 1e-4
par.Emax = 1.0
par.ρ0 = 1.0
par.tol = 1e-2
par.Neumann_bc = Ferrite.getfacetset(grid, "load")  # Nodes on the edge

γ = [1, 2, 3]                     # Penalty factor
vf =[0.0, 0.25, 0.5, 0.75]                   # Volume fraction
η = [π/(3.2) , π /(3.5) , π/(4.0)]               # Filter parameter
k = [4 , 8 , 12]                     # Sensitivity parameter


file_name = "optim"
dir = "/Users/aminalibakhshi/Desktop/pump_vtu/3d_cases/cantilever_3d/"
remove_vtk_files(dir)
par.max_itr = 200

@time optim_3D_combine(par, E, γ, vf, η, k, file_name, dir)