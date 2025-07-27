# Example for 2D Cantilever

Load the package
```
using Revise
using Ferrite
using PUPM
```
Generate a dynamic struc for fem inputs

```
par = DynamicParams() 
```
Function to generate grid (mesh and geometry)

```
# faces and nodes at the boundary are defiend
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "traction", x -> x[1] ≈ Lx && norm(x[2] - 0.5) <= 0.05)
    return grid
end

```

Function to create CellValues and FacetValues

```
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end

```

Function to create DofHandler (degrees of freedom )

```
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.close!(dh)
    return dh
end

```

Function to create Dirichlet boundary conditions

```
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

```

Define parameters for fem input

```
Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 120, 60   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

par.tnele = length(grid.cells)  # Total number of elements
par.grid = grid
# Create DOF handler and constraints
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh )

# Create CellValues and FacetValues
par.cell_values, par.facet_values = create_values()

# Define loads
pressure_value = 1.  # Example pressure in Pascals
par.loads = [LoadCondition("traction_load", [0.0, -pressure_value])]  # Load applied to the "traction" facet
# Neumann BC facet set
par.Neumann_bc = getfacetset(grid, "traction")
```

Define top input parameters

```
# Material properties
par.E0 = 1.0                
par.E = fill(par.E0, Ferrite.getncells(grid))  # Initialize Young's modulus for all cells
par.ν = 0.3 # Poisson's ratio 
# Optimization parameters
par.Emin = 1e-4             # Minimum Young's moduluspar.
par.Emax = 1.0              # Maximum Young's modulus
par.ρ0 = 1.0                # Initial density
par.tol = 1e-3            # Convergence tolerance
par.γ = 1           # Penalty factor
par.η = π /(4)              # Filter parameter
par.k = 8                  # Sensitivity parameter
par.vf = 0.5                # Volume fraction
par.max_itr = 200
```
Define file name to save data in vtu and defin dir

```
file_name = "linear_elasticty"
dir = "/Users/aminalibakhshi/Desktop/vtu_geo/" #(change to your dir)
``` 

Function to remove already present vtu files in dir (optional)

```
remove_vtk_files(dir)

```
Function to run topology optimization
```
top_upm(par, file_name, dir)
```

## Plain program

```
using Revise
using Ferrite
using PUPM

par = DynamicParams() 

# faces and nodes at the boundary are defiend
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "traction", x -> x[1] ≈ Lx && norm(x[2] - 0.5) <= 0.05)
    return grid
end

function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end

function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.close!(dh)
    return dh
end

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 120, 60   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

par.tnele = length(grid.cells)  # Total number of elements
par.grid = grid
# Create DOF handler and constraints
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh )

# Create CellValues and FacetValues
par.cell_values, par.facet_values = create_values()

# Define loads
pressure_value = 1.  # Example pressure in Pascals
par.loads = [LoadCondition("traction_load", [0.0, -pressure_value])]  # Load applied to the "traction" facet
# Neumann BC facet set
par.Neumann_bc = getfacetset(grid, "traction")

# Material properties
par.E0 = 1.0                
par.E = fill(par.E0, Ferrite.getncells(grid))  # Initialize Young's modulus for all cells
par.ν = 0.3 # Poisson's ratio 
# Optimization parameters
par.Emin = 1e-4             # Minimum Young's moduluspar.
par.Emax = 1.0              # Maximum Young's modulus
par.ρ0 = 1.0                # Initial density
par.tol = 1e-3            # Convergence tolerance
par.γ = 1           # Penalty factor
par.η = π /(4)              # Filter parameter
par.k = 8                  # Sensitivity parameter
par.vf = 0.5                # Volume fraction
par.max_itr = 200

file_name = "linear_elasticty"
dir = "/Users/aminalibakhshi/Desktop/vtu_geo/" # (change to your dir)

remove_vtk_files(dir)

top_upm(par, file_name, dir)

```
