using LinearAlgebra
using CSV
using DataFrames
using Random
using Dates

include("generate_matrices2D.jl")

# Edit these defaults, or override them from the command line with key=value:
# julia PTB_chern_BHZ_new_random_joint_LDOS_disorder_chern.jl L=40 proportion_sites=0.33 site_seed=1001 disorder_seed=2101
t0 = 1.0
t = 1.0
m_0 = 1.0
x_periodic = 1
y_periodic = 1
L = 40
proportion_sites = 0.33
disorder_strength = 0.1

output_root = "data/chern_BHZ_disorder"
seed = nothing
site_seed = nothing
disorder_seed = nothing
run_id = nothing

function parse_command_line_overrides!()
    global t0, t, m_0, x_periodic, y_periodic, L, proportion_sites, disorder_strength, output_root, seed, run_id
    global site_seed, disorder_seed

    for arg in ARGS
        if !occursin("=", arg)
            error("Arguments must be key=value, got: $arg")
        end

        key, value = split(arg, "=", limit=2)
        if key == "t0"
            t0 = parse(Float64, value)
        elseif key == "t"
            t = parse(Float64, value)
        elseif key == "m_0"
            m_0 = parse(Float64, value)
        elseif key == "x_periodic"
            x_periodic = parse(Int64, value)
        elseif key == "y_periodic"
            y_periodic = parse(Int64, value)
        elseif key == "L"
            L = parse(Int64, value)
        elseif key == "proportion_sites"
            proportion_sites = parse(Float64, value)
        elseif key == "disorder_strength"
            disorder_strength = parse(Float64, value)
        elseif key == "output_root"
            output_root = value
        elseif key == "seed"
            seed = parse(Int64, value)
        elseif key == "site_seed"
            site_seed = parse(Int64, value)
        elseif key == "disorder_seed"
            disorder_seed = parse(Int64, value)
        elseif key == "run_id"
            run_id = value
        else
            error("Unknown parameter: $key")
        end
    end
end

function validate_inputs(L, x_periodic, y_periodic, proportion_sites)
    L > 2 || error("L must be larger than 2 so the random PTB has interior sites.")
    x_periodic in (0, 1) || error("x_periodic must be 0 or 1.")
    y_periodic in (0, 1) || error("y_periodic must be 0 or 1.")
    0 < proportion_sites <= 1 || error("proportion_sites must be in the interval (0, 1].")
end

function lattice_points(Lx, Ly)
    points_x_array = zeros(Lx * Ly)
    points_y_array = zeros(Lx * Ly)

    for jj in 1:Ly
        for ii in 1:Lx
            index = (jj - 1) * Lx + ii
            points_x_array[index] = ii
            points_y_array[index] = jj
        end
    end

    return points_x_array, points_y_array
end

function find_random_non_nn_sites(rng, N_target, Lx, Ly, points_x_array, points_y_array)
    PTB_index_bigger = Int[]
    for jj in 1:Ly
        for ii in 1:Lx
            if (ii + jj) % 2 == 1
                push!(PTB_index_bigger, (jj - 1) * Lx + ii)
            end
        end
    end

    if N_target > length(PTB_index_bigger)
        error("Requested $N_target PTB sites, but only $(length(PTB_index_bigger)) checkerboard sites are available.")
    end

    PTB_index_bigger = shuffle(rng, shuffle(rng, PTB_index_bigger))
    return PTB_index_bigger[1:N_target]
end

"""Return the site coordinates and summed LDOS of the two states closest to zero energy."""
function generate_joint_LDOS(
    energy_eigenvalues_PTB,
    eigenstates_PTB,
    PTB_index,
    points_x_array,
    points_y_array,
)
    N_PTB = length(PTB_index)
    size(eigenstates_PTB, 1) == 2 * N_PTB || error("The eigenstate and PTB-site dimensions do not agree.")

    closest_to_zero = sortperm(abs.(energy_eigenvalues_PTB))[1:2]
    joint_LDOS_array = zeros(Float64, N_PTB)

    for ii in 1:N_PTB
        for state_index in closest_to_zero
            joint_LDOS_array[ii] += abs2(eigenstates_PTB[2 * ii - 1, state_index])
            joint_LDOS_array[ii] += abs2(eigenstates_PTB[2 * ii, state_index])
        end
    end

    points_PTB_array_x = points_x_array[PTB_index]
    points_PTB_array_y = points_y_array[PTB_index]

    return DataFrame(
        points_PTB_array_x=points_PTB_array_x,
        points_PTB_array_y=points_PTB_array_y,
        joint_LDOS_array=joint_LDOS_array,
    )
end

function generate_local_chern_marker(P_PTB, Q_PTB, XList_PTB, YList_PTB, PTB_index)
    N_PTB = length(PTB_index)
    size(P_PTB) == (2 * N_PTB, 2 * N_PTB) || error("The projector and PTB-site dimensions do not agree.")

    X1_PTB = Diagonal(kron(XList_PTB, [1, 1]))
    Y1_PTB = Diagonal(kron(YList_PTB, [1, 1]))

    # Bianco-Resta local Chern marker, as implemented in
    # PTB_chern_BHZ_new_random.ipynb.
    chern_matrix_PTB = -4 * pi * imag(P_PTB * X1_PTB * Q_PTB * Y1_PTB * P_PTB)
    chern_orbital_diagonal = diag(chern_matrix_PTB)
    local_chern_marker = [
        chern_orbital_diagonal[2 * ii - 1] + chern_orbital_diagonal[2 * ii]
        for ii in 1:N_PTB
    ]

    return DataFrame(
        site_index=1:N_PTB,
        lattice_index=PTB_index,
        x=XList_PTB,
        y=YList_PTB,
        local_chern_marker=local_chern_marker,
    )
end

function calculate_bott_index(t0, t, m_0, x_periodic, y_periodic, L, proportion_sites, disorder_strength;
    site_rng=Random.default_rng(), disorder_rng=Random.default_rng())
    Lx = L
    Ly = L
    N_PTB_requested = Int(round(Lx * Ly * proportion_sites))

    points_x_array, points_y_array = lattice_points(Lx, Ly)
    PTB_index = find_random_non_nn_sites(
        site_rng,
        N_PTB_requested,
        Lx,
        Ly,
        points_x_array,
        points_y_array)

    N_PTB = length(PTB_index)
    N_PTB > 0 || error("No PTB sites were selected. Increase L or proportion_sites.")

    println("Number of sites in PTB = ", N_PTB)
    println("Amount of sites in PTB = ", 100 * N_PTB / L^2, " %")
    println("PTB_index = ", PTB_index)

    sigma_x = [0 1; 1 0]
    sigma_y = [0 -im; im 0]
    sigma_z = [1 0; 0 -1]

    const_2D, CX2D, SX2D, CY2D, SY2D, _, _, _ = generate_matrices_2D(Lx, Ly, x_periodic, y_periodic)

    h_SC_chern = t * kron(SX2D, sigma_x) +
                 t * kron(SY2D, sigma_y) +
                 kron(m_0 * const_2D - t0 * (CX2D + CY2D), sigma_z)

    disorder_values = disorder_strength .* (2 .* rand(disorder_rng, Lx * Ly) .- 1)
    disorder_values .-= sum(disorder_values) / length(disorder_values)

    for ii in 1:(Lx * Ly)
        h_SC_chern[2 * ii - 1, 2 * ii - 1] += disorder_values[ii]
        h_SC_chern[2 * ii, 2 * ii] += disorder_values[ii]
    end

    if Hermitian_Check(h_SC_chern) == false
        error("Hamiltonian is not Hermitian")
    end

    PTB_orbital_index = Vector{Int64}(undef, 2 * N_PTB)
    for ii in 1:N_PTB
        PTB_orbital_index[2 * ii - 1] = 2 * PTB_index[ii] - 1
        PTB_orbital_index[2 * ii] = 2 * PTB_index[ii]
    end

    outside_orbital_index = setdiff(1:(2 * Lx * Ly), PTB_orbital_index)
    NOrbitalsInside = 2 * N_PTB
    NOrbitalsOutside = 2 * Lx * Ly - NOrbitalsInside

    H_11 = h_SC_chern[PTB_orbital_index, PTB_orbital_index]
    H_22 = h_SC_chern[outside_orbital_index, outside_orbital_index]
    H_12 = h_SC_chern[PTB_orbital_index, outside_orbital_index]
    H_21 = h_SC_chern[outside_orbital_index, PTB_orbital_index]

    H_22_regularized = Hermitian(H_22 + 1e-8 * Matrix(1.0I, NOrbitalsOutside, NOrbitalsOutside))
    H_PTB_renor = H_11 - Hermitian(H_12 * (H_22_regularized \ H_21))

    if Hermitian_Check(H_PTB_renor) == false
        error("Renormalized PTB Hamiltonian is not Hermitian")
    end

    YList_PTB = [cld(index, Lx) for index in PTB_index]
    XList_PTB = PTB_index .- (YList_PTB .- 1) .* Lx

    XListKron_PTB = kron(XList_PTB, [1, 1])
    YListKron_PTB = kron(YList_PTB, [1, 1])

    U_X = diagm(exp.(im * 2 * pi .* XListKron_PTB ./ Lx))
    V_Y = diagm(exp.(im * 2 * pi .* YListKron_PTB ./ Ly))

    energy_eigenvalues_PTB, eigenstates_PTB = eigen(H_PTB_renor)
    joint_LDOS_data = generate_joint_LDOS(
        energy_eigenvalues_PTB,
        eigenstates_PTB,
        PTB_index,
        points_x_array,
        points_y_array,
    )
    filled_eigenstates_PTB = eigenstates_PTB[:, 1:N_PTB]

    P_PTB = conj(filled_eigenstates_PTB) * transpose(filled_eigenstates_PTB)
    Q_PTB = Matrix(1.0I, 2 * N_PTB, 2 * N_PTB) - P_PTB

    local_chern_data = generate_local_chern_marker(
        P_PTB, Q_PTB, XList_PTB, YList_PTB, PTB_index,
    )

    U = P_PTB * U_X * P_PTB + Q_PTB
    V = P_PTB * V_Y * P_PTB + Q_PTB

    bott = U * V * U' * V'
    eigvals_bott = eigvals(bott)
    bott_index = real((-im / (2 * pi)) * sum(log.(Complex.(eigvals_bott))))

    gap_PTB = 2 * minimum(abs.(energy_eigenvalues_PTB))
    println("Gap_PTB = ", gap_PTB)
    println("Bott index = ", bott_index)
    println("size H_PTB_renor = ", size(H_PTB_renor))

    bott_index_data = DataFrame(
        t0=[t0],
        t=[t],
        m_0=[m_0],
        x_periodic=[x_periodic],
        y_periodic=[y_periodic],
        L=[L],
        proportion_sites=[proportion_sites],
        disorder_strength=[disorder_strength],
        n_ptb=[N_PTB],
        gap_PTB=[gap_PTB],
        bott_index=[bott_index],
    )

    return bott_index_data, joint_LDOS_data, local_chern_data
end

function output_filename(t0, t, m_0, x_periodic, y_periodic, L, proportion_sites, disorder_strength, run_id;
    seed=nothing, site_seed=nothing, disorder_seed=nothing)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    stem = "t0=$(t0)_t=$(t)_m_0=$(m_0)_x_periodic=$(x_periodic)_y_periodic=$(y_periodic)_L=$(L)_proportion_sites=$(proportion_sites)_disorder_strength=$(disorder_strength)"
    for (name, value) in (("seed", seed), ("site_seed", site_seed), ("disorder_seed", disorder_seed))
        if value !== nothing
            stem *= "_$(name)=$(value)"
        end
    end
    if run_id !== nothing
        stem = "$(stem)_run_id=$(run_id)"
    end
    return "$(stem)_$(timestamp).csv"
end

function main()
    parse_command_line_overrides!()
    validate_inputs(L, x_periodic, y_periodic, proportion_sites)

    for (name, value) in (("seed", seed), ("site_seed", site_seed), ("disorder_seed", disorder_seed))
        value === nothing || value >= 0 || error("$name must be nonnegative.")
    end
    # Preserve legacy seed= behavior; explicit seeds override each stream independently.
    if seed !== nothing
        Random.seed!(seed)
    end
    site_rng = site_seed === nothing ? Random.default_rng() : MersenneTwister(site_seed)
    disorder_rng = disorder_seed === nothing ? Random.default_rng() : MersenneTwister(disorder_seed)

    t_start = time()
    println("t0 = ", t0)
    println("t = ", t)
    println("m_0 = ", m_0)
    println("x_periodic = ", x_periodic)
    println("y_periodic = ", y_periodic)
    println("L = ", L)
    println("proportion_sites = ", proportion_sites)
    println("disorder_strength = ", disorder_strength)
    println("seed = ", seed)
    println("site_seed = ", site_seed)
    println("disorder_seed = ", disorder_seed)

    bott_index_data, joint_LDOS_data, local_chern_data = calculate_bott_index(
        t0,
        t,
        m_0,
        x_periodic,
        y_periodic,
        L,
        proportion_sites,
        disorder_strength;
        site_rng=site_rng,
        disorder_rng=disorder_rng,
    )

    for data in (bott_index_data, joint_LDOS_data, local_chern_data)
        data[!, :seed] = fill(something(seed, missing), nrow(data))
        data[!, :site_seed] = fill(something(site_seed, missing), nrow(data))
        data[!, :disorder_seed] = fill(something(disorder_seed, missing), nrow(data))
    end
    filename = output_filename(t0, t, m_0, x_periodic, y_periodic, L, proportion_sites, disorder_strength, run_id;
        seed=seed, site_seed=site_seed, disorder_seed=disorder_seed)

    output_dir = joinpath(output_root, "bott_index")
    mkpath(output_dir)

    path = joinpath(output_dir, filename)
    CSV.write(path, bott_index_data)

    println("Exported bott_index to ", path)

    joint_LDOS_output_dir = joinpath(output_root, "joint_LDOS")
    mkpath(joint_LDOS_output_dir)
    joint_LDOS_path = joinpath(
        joint_LDOS_output_dir,
        filename,
    )
    CSV.write(joint_LDOS_path, joint_LDOS_data)

    println("Exported joint LDOS to ", joint_LDOS_path)
    local_chern_output_dir = joinpath(output_root, "chern_marker_sitewise")
    mkpath(local_chern_output_dir)
    local_chern_path = joinpath(local_chern_output_dir, filename)
    CSV.write(local_chern_path, local_chern_data)
    println("Exported local Chern marker to ", local_chern_path)
    println("time taken = ", time() - t_start, " seconds")
end

main()
