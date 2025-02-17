using Revise
using Ferrite
using PUPM
### 3d topology optimization for chair
@time begin
par = DynamicParams()
# Function to check if a point is inside a circle on a given plane
function in_circle(x, r, cx, cy, cz)
    return (x[3] ≈ cz) && ((x[1] - cx)^2 + (x[2] - cy)^2 <= r^2)
end

# Function to create a 3D grid
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    # Define the domain boundaries using Vec{3}
    left = Vec(0.0, 0.0, 0.0)
    right = Vec(Lx, Ly, Lz)

    # Generate the 3D hexahedral grid
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    return grid
end

# Function to define boundary conditions and nodesets
function create_boundary(grid, Lx, Ly, Lz)
    # Parameters for the circle on the top surface
    r_top = 20.0
    cx_top, cy_top = Lx / 2, Ly / 2

    # Parameters for the circles on the bottom surface
    r_bottom = 10.0

    # Define the bottom corner centers
    corner_centers = [
        (0.0, 0.0),
        (Lx, 0.0),
        (0.0, Ly),
        (Lx, Ly)
    ]

    # Add nodeset for the top surface inside the circle
    # addnodeset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, Lz))
    addfacetset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, Lz))

    # Add nodesets for each bottom corner circle
    for (i, (cx, cy)) in enumerate(corner_centers)
        addnodeset!(grid, "bottom_corner_circle_$i", x -> in_circle(x, r_bottom, cx, cy, 0.0))
    end
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
    ch = Ferrite.ConstraintHandler(dh)
    for i in 1:4
        dbc = Dirichlet(:u, getnodeset(grid, "bottom_corner_circle_$i"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
        add!(ch, dbc)
    end
    Ferrite.close!(ch)
    return ch
end


# Main script
Lx, Ly, Lz = 100.0, 100.0, 100.0
nx, ny, nz = 15, 15, 15

grid = create_grid(Lx, Ly, Lz, nx, ny, nz)
par.tnele = length(grid.cells)  # Total number of elements
create_boundary(grid, Lx, Ly, Lz)
par.grid = grid
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh, grid)

par.cell_values, par.facet_values = create_values()


#par.loads = [LoadCondition_3d("nodal_load", [0.0, 0.0, 1.0])]
par.loads = [LoadCondition_3d("traction_load", [0.0, 0.0, 1.0])]


# Material properties
par.E0 = 1.0
par.E = fill(par.E0,Ferrite.getncells(grid))
par.ν = 0.3

# Optimization parameters
par.Emin = 1e-4
par.Emax = 1.0
par.ρ0 = 1.0
par.tol = 1e-3
par.γ = 1
par.η = π / (3.0)
par.k = 2
par.vf = 0.5

# Neumann BC
# par.Neumann_bc = Ferrite.getnodeset(grid, "top_circle")  # Nodes on the edge
par.Neumann_bc = Ferrite.getfacetset(grid, "top_circle")  # Nodes on the edge


file_name = "linear_elasticity_3d"
dir = "/Users/aminalibakhshi/Desktop/pump_vtu"
par.max_itr = 300
remove_vtk_files(dir)
# fem example
#fem = fem_solver_3d(par) 
# Run topology optimization
top_upm_3d(par, file_name, dir)

end 