using Revise
using Ferrite
using PUPM
@time begin
par = DynamicParams() 
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0) 
    addnodeset!(grid, "support_2", x -> x[2] ≈ 0.0)  
    addfacetset!(grid, "pressure", x -> x[1] ≈ Lx) 
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
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0], [1]))
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_2"), (x, t) -> [0.0], [2]))
    Ferrite.close!(ch)
    return ch
end
Lx, Ly = 1.0, 1.0  
nx, ny = 10, 10   
par.grid = create_grid(Lx, Ly, nx, ny)
par.dh = create_dofhandler(par.grid)
par.ch = create_bc(par.dh)
par.cell_values, par.facet_values = create_values()
pressure_value = 1e10  # Example pressure in Pascals
par.loads = [LoadCondition("pressure_load", pressure_value)] 
par.E = fill(210e9, nx * ny)
par.ν = 0.3 
dh = par.dh; grid = dh.grid
par.Neumann_bc = Ferrite.getfacetset(dh.grid, "pressure")
result = fem_solver(par)
u = result.u
maximum(u)
end


