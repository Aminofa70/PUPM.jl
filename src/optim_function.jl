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
function filter_density_to_vf!(density, vf, tnele, eta)
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
            sumdmin += transform(density[i], rhominbound, eta, rhomin, rhomax) / (tnele)
        end
        sumdmax = 0.0
        for i in eachindex(density)
            sumdmax += transform(density[i], rhomaxbound, eta, rhomin, rhomax) / (tnele)
        end
        sumdmid = 0.0
        for i in eachindex(density)
            sumdmid += transform(density[i], rhotr, eta, rhomin, rhomax) / (tnele)
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
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
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
        ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
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
        E_node = fem.E_node
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
            write_node_data(vtk, E_node, "Nodal Young's modulus")
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
    E_node = fem.E_node
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
        write_node_data(vtk, E_node, "Nodal Young's modulus")
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
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
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
        ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
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
        E_node = fem.E_node
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
            write_node_data(vtk, E_node, "Nodal Young's modulus")
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
    E_node = fem.E_node
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
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end
end