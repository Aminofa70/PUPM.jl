using Revise
using Ferrite
using PUPM
### topology optmization MBB beam
# Function to create the grid
@time begin
par = DynamicParams() 
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0) #fixed in x-direction
    addnodeset!(grid, "support_2", x -> x[1] ≈ Lx && x[2] ≈ 0.0) # fixed in y direction
    addnodeset!(grid, "nodal_force", x -> x[1] ≈ 0.0 && x[2] ≈ Ly) # nodal_force
    return grid
end

# Function to create CellValues and FacetValues
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create DofHandler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.close!(dh)
    return dh
end

# Function to create Dirichlet boundary conditions
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "support_1"), (x, t) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "support_2"), (x, t) -> [0.0], [2]))
    Ferrite.close!(ch)
    return ch
end

# Define parameters for the plate and mesh
Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 120, 60   # Number of elements along x and y

grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
par.tnele = length(grid.cells)  # Total number of elements

par.grid = grid
# Create DOF handler and constraints
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh)

# Create CellValues and FacetValues
par.cell_values, par.facet_values = create_values()

# Define loads
par.loads = [LoadCondition("nodal_load", [0.0, +1.0])]  # Load applied to the "traction" facet

# Material properties
par.E0 = 1.0                # Initial Young's modulus
par.E = fill(par.E0, Ferrite.getncells(grid))  # Initialize Young's modulus for all cells
par.ν = 0.3  #  Poisson's ratio 
# Optimization parameters
par.Emin = 1e-4             # Minimum Young's modulus
par.Emax = 1.0              # Maximum Young's modulus
par.ρ0 = 1.0                # Initial density
par.tol = 1e-3            # Convergence tolerance
par.γ = 3               # Penalty factor
par.η = π / (3.5)              # Filter parameter
par.k = 8                   # Sensitivity parameter
par.vf = 0.5          # Volume fraction
par.max_itr = 200
# Neumann BC 
par.Neumann_bc = "nodal_force"
file_name = "linear_elasticty"
dir = "/Users/aminalibakhshi/Desktop/data_vtu"
remove_vtk_files(dir)
# Run the topology optimization
top_upm!(par, file_name, dir)


end 
