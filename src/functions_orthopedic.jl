"""
stiffness matrix for orthotropic

"""
function C_orthotropic(Ex, Ey, νxy)

    νyx = (Ey * νxy) / Ex
    Q11 = Ex / (1 - νxy * νyx)
    Q12 = (νxy * Ey) / (1 - νxy * νyx)

    Q21 = (νxy * Ey) / (1 - νxy * νyx)
    Q22 = Ey / (1 - νxy * νyx)

    Gxy = (Ex * Ey) / (Ex + Ey + 2 * Ey * νxy)
    Q66 = Gxy

    C_voigt = [Q11 Q12 0.0;
        Q21 Q22 0.0;
        0.0 0.0 Q66]

    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end # end of function
####################################################
function compute_derivatives(Ex, Ey, νxy)
    # Define a function that takes a vector of parameters [Ex, Ey] and returns the flattened C_orthotropic tensor
    function f(params)
        C = C_orthotropic(params[1], params[2], νxy)
        return tovoigt(C)  # Convert the tensor to Voigt notation (3x3 matrix)
    end

    # Compute the Jacobian using ForwardDiff
    params = [Ex, Ey]
    jacobian = ForwardDiff.jacobian(f, params)

    # Reshape the Jacobian to match the tensor structure and convert to SymmetricTensor{4,2}
    dC_dEx = fromvoigt(SymmetricTensor{4,2}, reshape(jacobian[:, 1], (3, 3)))  # Derivative with respect to Ex
    dC_dEy = fromvoigt(SymmetricTensor{4,2}, reshape(jacobian[:, 2], (3, 3)))  # Derivative with respect to Ey

    return dC_dEx, dC_dEy
end

####################################################
function assemble_cell_orthotropic!(ke, cell_values, Ex, Ey, νxy)
    C = C_orthotropic(Ex, Ey, νxy)
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
####################################################
function assemble_global_orthotropic!(K, dh, cell_values, Ex, Ey, νxy)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)

        assemble_cell_orthotropic!(ke, cell_values, Ex[cell_index], Ey[cell_index], νxy)
        assemble!(assembler, celldofs(cell), ke)
    end

    return K
end
####################################################
function calculate_strain_energy_orthotropic(grid, dh, cv, u, Ex, Ey, νxy)
    element_strain_energies = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)
        
        C = C_orthotropic(Ex[cell_index], Ey[cell_index], νxy)
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

#######################################################
function dC_orthotropic_dEx(Ex, Ey, νxy)
    x = Ex
    y = Ey
    r = νxy
    dQ11_dEx = (-r^2 * y) / (x * (1 - (r^2 * y) / x)^2) + (1 / (1 - (r^2 * y) / x))
    dQ12_dEx = (-r^3 * y^2) / (x^2 * (1 - (r^2 * y) / x)^2)
    dQ21_dEx = (-r^3 * y^2) / (x^2 * (1 - (r^2 * y) / x)^2)
    dQ22_dEx = (-r^2 * y^2) / (x^2 * (1 - (r^2 * y) / x)^2)

    dQ66_dEx = (-x * y) / (x + y + 2 * r * y)^2 + (y / (x + y + 2 * r * y))
    dC_voigt = [dQ11_dEx dQ12_dEx 0.0;
                dQ21_dEx dQ22_dEx 0.0;
                0.0 0.0 dQ66_dEx]
    
    return fromvoigt(SymmetricTensor{4,2}, dC_voigt)
end

function dC_orthotropic_dEy(Ex, Ey, νxy)
    
    x = Ex
    y = Ey
    r = νxy
    dQ11_dEy = r^2 / (1 - (r^2 * y) / x)^2
    dQ12_dEy = (r^3 * y) / (x * (1 - (r^2 * y) / x)^2) + (r / (1 - (r^2 * y) / x))
    dQ21_dEy = (r^3 * y) / (x * (1 - (r^2 * y) / x)^2) + (r / (1 - (r^2 * y) / x))
    dQ22_dEy = (r^2 * y) / (x * (1 - (r^2 * y) / x)^2) + (1 / (1 - (r^2 * y) / x))

    dQ66_dEy = (-(1 + 2 * r) * x * y) / (x + y + 2 * r * y)^2 + (x / (x + y + 2 * r * y)) 
   

    dC_voigt = [dQ11_dEy dQ12_dEy 0.0;
                dQ21_dEy dQ22_dEy 0.0;
                0.0 0.0 dQ66_dEy]
    
    return fromvoigt(SymmetricTensor{4,2}, dC_voigt)
end



##################################################
function calculate_Hx(grid, dh, cv, u, Ex, Ey, νxy)
    element_strain_energy_derivatives = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)
        
        # Handle E as either scalar or vector
        dC_dE = dC_orthotropic_dEx(Ex[cell_index], Ey[cell_index], νxy)
        #dC_dE, _ = compute_derivatives(Ex[cell_index], Ey[cell_index], νxy)
        
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
function calculate_Hy(grid, dh, cv, u, Ex, Ey, νxy)
    element_strain_energy_derivatives = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cv, cell)

        # Handle E as either scalar or vector
        dC_dE = dC_orthotropic_dEy(Ex[cell_index], Ey[cell_index], νxy)
        # _, dC_dE = compute_derivatives(Ex[cell_index], Ey[cell_index], νxy)
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


struct FEMSolver_orthotropic
    u::Vector{Float64}
    compliance::Float64
    U::Vector{Float64}
    Hx::Vector{Float64}
    Hy::Vector{Float64}
end

function fem_solver_combine_orthotropic(par::DynamicParams, Ex, Ey)
    # Extract common parameters from the input
    grid = par.grid
    cell_values = par.cell_values
    facet_values = par.facet_values
    dh = par.dh
    ch = par.ch
    # E = par.E
    νxy = par.νxy

    # Allocate and assemble the global stiffness matrix
    K = allocate_matrix(dh)
    assemble_global_orthotropic!(K, dh, cell_values, Ex, Ey, νxy)

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
    
    U = calculate_strain_energy_orthotropic(grid, dh, cell_values, u, Ex, Ey, νxy)
    Hx = calculate_Hx(grid, dh, cell_values, u, Ex, Ey, νxy)
    Hy = calculate_Hy(grid, dh, cell_values, u, Ex, Ey, νxy)

    # Return the result
    return FEMSolver_orthotropic(u, c, U, Hx,Hy)
end


function compute_gamma(E::Vector, E0::Real, ρ::Vector, ρ0::Real, γ::Real)
    # Check if E and ρ have the same length
    if length(E) != length(ρ)
        error("Vectors E and ρ must have the same length.")
    end

    # Compute γ_y for each pair (E, ρ)
    γ_y = log.(E ./ E0) ./ log.(ρ ./ ρ0)

    # Clamp γ_y to a specified range (e.g., [1.0, 5.0])
    γ_min = 1.0
    γ_max = 5.0
    γ_y = clamp.(γ_y, γ_min, γ_max)

    # Replace NaN values in γ_y with the fallback scalar γ
    γ_y[isnan.(γ_y)] .= γ

    return γ_y
end

function transfer_to_young_y_dir(ρnew, E0,ρ0, γ, Emin, Emax)

    # Ensure ρnew and γ have the same length
    @assert length(ρnew) == length(γ) "ρnew and γ must have the same length"

    # Calculate Enew element-wise
    Enew = E0 * (ρnew / ρ0) .^ γ
    Enew = [clamp(x, Emin, Emax) for x in Enew]

    return Enew
end
###############################################################
function top_upm_orthotropic(par::DynamicParams, Ex, Ey,  k, γ, η ,volfrac, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    #E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    #k = par.k ; γ = par.γ ; volfrac = par.vf; 
    #η = par.η; 
    ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    ## Initial FEM solve
    fem = fem_solver_combine_orthotropic(par, Ex, Ey)
    compliance = fem.compliance

    Hx = fem.Hx
    Hy = fem.Hy
    W_tot = sum(fem.U)
    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

    ## println("Iter $loop: C = $compliance, ΔR = $A")

    ## Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        fem = fem_solver_combine_orthotropic(par, Ex, Ey)
        compliance = fem.compliance
        Hx = fem.Hx
        Hy = fem.Hy
        W_tot = sum(fem.U)
        
        if volfrac == 0.0
            
            Enew_x = update_upm(k, Ex, Hx, Emax, Emin)
            Enew_frac_x = Enew_x

            Enew_y = update_upm(k, Ey, Hy, Emax, Emin)
            Enew_frac_y = Enew_y

        elseif volfrac > 0.0

            # Enew_x = update_upm(k, Ex, Hx, Emax, Emin)
            # ρ = transfer_to_density(Enew_x, E0, ρ0, γ)
            # ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            # Enew_frac_x = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
            # #-----------------------------------
            # Enew_y = update_upm(k, Ey, Hy, Emax, Emin)
            # Enew_frac_y = Enew_y

            Enew_x = update_upm(k, Ex, Hx, Emax, Emin)
            ρ = transfer_to_density(Enew_x, E0, ρ0, γ)
            Enew_y = update_upm(k, Ey, Hy, Emax, Emin)
            γ_y = compute_gamma(Enew_y, E0, ρ, ρ0, γ)
            ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
            Enew_frac_x = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)
            Enew_y = transfer_to_young_y_dir(ρnew, E0, ρ0, γ_y, Emin, Emax)
            Enew_frac_y = Enew_y

        else
            error("Invalid value for volfrac")
        end

        # Update E in par so that fem_solver uses the updated material distribution
        # par.E = Enew_frac
        Ex = Enew_frac_x
        Ey = Enew_frac_y

        fem = fem_solver_combine_orthotropic(par, Ex, Ey)
        compliance = fem.compliance
        Hx = fem.Hx
        Hy = fem.Hy
        W_tot = sum(fem.U)
       
        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        ## println("Iter $loop: C = $compliance, ΔR = $A")

    end

    ## Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final write
    fem = fem_solver_combine_orthotropic(par, Ex, Ey)
    compliance = fem.compliance
    u = fem.u
    U = fem.U

    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        
        write_cell_data(vtk, Ex, "Young's modulus Ex")
        write_cell_data(vtk, Ey, "Young's modulus Ey")
        write_cell_data(vtk, U, "Strain Energy")
        Ferrite.write_cellset(vtk, grid)
    end

    return compliance
end

function optim_2D_combine_orthotropic(
    par::DynamicParams, 
    Ex,
    Ey, 
    gamma_values, 
    volfrac_values, 
    eta_values, 
    k_values, 
    name_of_file::String, 
    directory::String
)
    # Ensure output directory exists
    if !isdir(directory)
        println("Directory does not exist. Creating: $directory")
        mkpath(directory)
    end

    # Initialize log file
    log_path = joinpath(directory, string(name_of_file, ".txt"))
    open(log_path, "w") do file
        write(file, "File γ vf η k Compliance\n")  # Space-separated header
    end

    index = 1  # Initialize index for file naming and logging

    for vf in volfrac_values
        for γ in gamma_values
            for k in k_values
                η_values = vf == 0.0 ? [0.0] : eta_values  # Single η value when vf == 0.0
                nvf_str = vf == 0.0 ? "nvf" : string(vf)   # Use "nvf" for zero volume fraction
                
                for η in η_values
                    file_name = string("out_", lpad(index, 4, '0'))  # Ensure unique filename
                    file_name_with_params = joinpath(directory, string(file_name, ".vtu"))

                    # Compute compliance
                    compliance = top_upm_orthotropic(par, Ex,Ey, k, γ, η, vf, file_name_with_params, directory)

                    # Append results to log file
                    open(log_path, "a") do file
                        write(file, "$file_name $γ $nvf_str $η $k $compliance\n")
                    end

                    println("Saved file: $file_name_with_params with compliance: $compliance")
                    index += 1
                end
            end
        end
    end

    println("Optimization completed. Results saved to $log_path")
end