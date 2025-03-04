using Revise
using Ferrite
using PUPM
# test fem_solver
"""
Finite Element programming for 2D plane stress 
Plate under tension
Reference:
https://shop.elsevier.com/books/practical-programming-of-finite-element-procedures-for-solids-and-structures-with-matlab/farahmand-tabar/978-0-443-15338-9

"""
par = DynamicParams()  # Create a dynamic parameter object
# Function to create the grid
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)  # left side; Dirichlet BC
    addnodeset!(grid, "support_2", x -> x[2] ≈ 0.0)  # bottom side; Dirichlet BC
    addfacetset!(grid, "pressure", x -> x[1] ≈ Lx)   # right side; Neumann BC
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
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0], [1]))
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_2"), (x, t) -> [0.0], [2]))
    Ferrite.close!(ch)
    return ch
end

# Define parameters
Lx, Ly = 1.0, 1.0  # Plate dimensions
nx, ny = 4, 4      # Number of elements along x and y
par.grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

# Create DOF handler and constraints
par.dh = create_dofhandler(par.grid)
par.ch = create_bc(par.dh)

# Create CellValues and FacetValues
par.cell_values, par.facet_values = create_values()

# # # Define loads
pressure_value = 1e10  # Example pressure in Pascals
par.loads = [LoadCondition("pressure_load", pressure_value)]  # Load applied to "pressure" facet

# Material properties
par.E = fill(210e9, nx * ny)  # Young's modulus (Pa)
#par.E = 210e9
par.ν = 0.3 # Poisson's ratio

# Neumann BC facet set
dh = par.dh; grid = dh.grid
par.Neumann_bc = Ferrite.getfacetset(dh.grid, "pressure")

# Solve the FEM problem using OptiUPM
result = fem_solver(par)

u = result.u

display(maximum(u))
display(minimum(u))