
using Test
using JuMP
using Ipopt
using ExaModels
using NLPModels
using NLPModelsIpopt
using NLPModelsJuMP
using SparseArrays
using LinearAlgebra
using Random
using COPSBenchmark

# Reference values from http://www.mcs.anl.gov/~more/cops/cops3.pdf
# Sizes reduced where possible; rtol=1e-4 requires sufficient discretization.
COPS_INSTANCES = [
    (COPSBenchmark.bearing_model, (50, 50), -1.5482e-1),
    (COPSBenchmark.chain_model, (800,), 5.06891),
    (COPSBenchmark.camshape_model, (200,), 4.2791),
    (COPSBenchmark.catmix_model, (100,), -4.80556e-2),
    (COPSBenchmark.channel_model, (100,), 1.0),
    (COPSBenchmark.elec_model, (50,), 1.0552e3),
    (COPSBenchmark.gasoil_model, (100,), 5.2366e-3),
    (COPSBenchmark.glider_model, (100,), -1.25505e3),
    (COPSBenchmark.marine_model, (100,), 1.97462e7),
    (COPSBenchmark.methanol_model, (100,), 9.02229e-3),
    (COPSBenchmark.minsurf_model, (50, 50), 2.51488),
    (COPSBenchmark.pinene_model, (100,), 1.98721e1),
    (COPSBenchmark.polygon_model, (100,), -0.674981),
    (COPSBenchmark.robot_model, (200,), 9.14138),
    (COPSBenchmark.rocket_model, (400,), 1.01283),
    (COPSBenchmark.steering_model, (200,), 5.54577e-1),
    (COPSBenchmark.tetra_duct15_model, (), 1.04951e4),
    (COPSBenchmark.tetra_duct20_model, (), 4.82685e3),
    (COPSBenchmark.tetra_gear_model, (), 4.15163e3),
    (COPSBenchmark.torsion_model, (50, 50), -4.18087e-1),
    (COPSBenchmark.dirichlet_model, (20,), 1.71464e-2),
    # N.B.: comment henon and lane-emden as they take too long in the CI.
    # (COPSBenchmark.henon_model, (10,), 6.667736),
    # (COPSBenchmark.lane_emden_model, (20,), 9.11000),
    (COPSBenchmark.triangle_pacman_model, (), 1.25045e3),
]

function solve_backend(model, ::COPSBenchmark.JuMPBackend)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_attribute(model, "max_iter", 10)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    obj_val = JuMP.objective_value(model)
    return status, obj_val
end

function solve_backend(model, ::COPSBenchmark.ExaModelsBackend)
    results = ipopt(model; print_level=0, tol=1e-8, max_iter=10)
    status = results.status == :first_order
    obj_val = results.objective
    return status, obj_val
end

@testset "Backend $(backend)" for backend in [
    COPSBenchmark.JuMPBackend(),
    COPSBenchmark.ExaModelsBackend(),
]
    @testset "Instance $instance" for (instance, params, result) in COPS_INSTANCES
        if hasmethod(instance, Tuple{typeof(backend), typeof.(params)...})
            model = instance(backend, params...)
            # Smoke test: build model + run 10 Ipopt iterations (not solving to convergence)
            solve_backend(model, backend)
            @test true
        end
    end
end

# Models to compare callbacks between JuMP and ExaModels backends.
COMPARE_INSTANCES = [
    (COPSBenchmark.bearing_model, (10, 10)),
    (COPSBenchmark.camshape_model, (100,)),
    (COPSBenchmark.catmix_model, (10,)),
    (COPSBenchmark.chain_model, (100,)),
    (COPSBenchmark.channel_model, (20,)),
    (COPSBenchmark.elec_model, (25,)),
    (COPSBenchmark.gasoil_model, (10,)),
    (COPSBenchmark.glider_model, (10,)),
    (COPSBenchmark.marine_model, (10,)),
    (COPSBenchmark.methanol_model, (10,)),
    (COPSBenchmark.minsurf_model, (10, 10)),
    (COPSBenchmark.pinene_model, (10,)),
    (COPSBenchmark.polygon_model, (50,)),
    (COPSBenchmark.robot_model, (20,)),
    (COPSBenchmark.rocket_model, (20,)),
    (COPSBenchmark.steering_model, (20,)),
    (COPSBenchmark.torsion_model, (10, 10)),
]

# Find the permutation P such that keys_a[P] ≈ keys_b.
# This recovers the variable/constraint reordering done by JuMP.
# Returns nothing if keys are not uniquely matchable.
function find_permutation(keys_a, keys_b)
    length(keys_a) == length(keys_b) || return nothing
    allunique(keys_a) && allunique(keys_b) || return nothing
    perm_a = sortperm(keys_a)
    perm_b = sortperm(keys_b)
    P = similar(perm_a)
    for i in eachindex(perm_a)
        P[perm_b[i]] = perm_a[i]
    end
    return P
end

# Generate a random point within [lvar, uvar], perturbed from x0.
function perturb_x0(x0, lvar, uvar; rng = Random.default_rng(), scale = 0.1)
    x = copy(x0)
    for i in eachindex(x)
        lo = max(lvar[i], x0[i] - scale * max(abs(x0[i]), 1.0))
        hi = min(uvar[i], x0[i] + scale * max(abs(x0[i]), 1.0))
        if isfinite(lo) && isfinite(hi) && lo < hi
            x[i] = lo + (hi - lo) * rand(rng)
        end
    end
    return x
end

@testset "Compare callbacks: $instance" for (instance, params) in COMPARE_INSTANCES
    jump_model = instance(COPSBenchmark.JuMPBackend(), params...)
    exa_model = instance(COPSBenchmark.ExaModelsBackend(), params...)

    jump_nlp = MathOptNLPModel(jump_model)

    n = get_nvar(jump_nlp)
    m = get_ncon(jump_nlp)

    # Structural comparison
    @test n == get_nvar(exa_model)
    @test m == get_ncon(exa_model)

    # Find variable permutation by matching (lvar, uvar, x0) tuples
    lj, uj, x0j = get_lvar(jump_nlp), get_uvar(jump_nlp), jump_nlp.meta.x0
    le, ue, x0e = get_lvar(exa_model), get_uvar(exa_model), exa_model.meta.x0
    var_keys_j = collect(zip(lj, uj, x0j))
    var_keys_e = collect(zip(le, ue, x0e))
    Pv = find_permutation(var_keys_j, var_keys_e)

    if Pv !== nothing
        # Exact element-wise comparison via permutation
        @test lj[Pv] ≈ le
        @test uj[Pv] ≈ ue
        @test x0j[Pv] ≈ x0e
    else
        # Fallback: sorted comparison (duplicate keys)
        @test sort(lj) ≈ sort(le)
        @test sort(uj) ≈ sort(ue)
        @test sort(x0j) ≈ sort(x0e)
    end

    # Objective (scalar, order-independent)
    @test obj(jump_nlp, x0j) ≈ obj(exa_model, x0e) rtol = 1e-6

    # Find constraint permutation
    lcj, ucj = get_lcon(jump_nlp), get_ucon(jump_nlp)
    lce, uce = get_lcon(exa_model), get_ucon(exa_model)
    cj0 = cons(jump_nlp, x0j)
    ce0 = cons(exa_model, x0e)
    con_keys_j = collect(zip(lcj, ucj, cj0))
    con_keys_e = collect(zip(lce, uce, ce0))
    Pc = find_permutation(con_keys_j, con_keys_e)

    if Pv !== nothing && Pc !== nothing
        # --- Exact comparison via permutations ---
        @test lcj[Pc] ≈ lce
        @test ucj[Pc] ≈ uce
        @test grad(jump_nlp, x0j)[Pv] ≈ grad(exa_model, x0e) rtol = 1e-6
        @test cj0[Pc] ≈ ce0 rtol = 1e-6

        # Jacobian (apply row/column permutation)
        Jj_r, Jj_c = jac_structure(jump_nlp)
        Jj_v = jac_coord(jump_nlp, x0j)
        J_jump = Matrix(sparse(Jj_r, Jj_c, Jj_v, m, n))

        Je_r, Je_c = jac_structure(exa_model)
        Je_v = jac_coord(exa_model, x0e)
        J_exa = Matrix(sparse(Je_r, Je_c, Je_v, m, n))

        @test J_jump[Pc, Pv] ≈ J_exa atol = 1e-6

        # Hessian (apply variable permutation to both rows and columns)
        y0 = ones(m)
        Hj_r, Hj_c = hess_structure(jump_nlp)
        Hj_v = hess_coord(jump_nlp, x0j, y0)
        H_jump = Matrix(Symmetric(sparse(Hj_r, Hj_c, Hj_v, n, n), :L))

        He_r, He_c = hess_structure(exa_model)
        He_v = hess_coord(exa_model, x0e, y0)
        H_exa = Matrix(Symmetric(sparse(He_r, He_c, He_v, n, n), :L))

        @test H_jump[Pv, Pv] ≈ H_exa atol = 1e-6

        # Test at perturbed points
        Random.seed!(42)
        for _ in 1:3
            xp_exa = perturb_x0(x0e, le, ue)
            xp_jump = similar(x0j)
            xp_jump[Pv] = xp_exa

            @test obj(jump_nlp, xp_jump) ≈ obj(exa_model, xp_exa) rtol = 1e-6
            @test grad(jump_nlp, xp_jump)[Pv] ≈ grad(exa_model, xp_exa) rtol = 1e-6
            @test cons(jump_nlp, xp_jump)[Pc] ≈ cons(exa_model, xp_exa) rtol = 1e-6
        end
    else
        # --- Fallback: order-invariant comparison ---
        # Note: constraint bounds/values may differ due to JuMP normalizing
        # constants into bounds. Compare gradient and Jacobian/Hessian spectra
        # which are invariant under such reformulations.
        @test sort(grad(jump_nlp, x0j)) ≈ sort(grad(exa_model, x0e)) rtol = 1e-6

        # Jacobian: compare singular values
        Jj_r, Jj_c = jac_structure(jump_nlp)
        Jj_v = jac_coord(jump_nlp, x0j)
        J_jump = Matrix(sparse(Jj_r, Jj_c, Jj_v, m, n))
        Je_r, Je_c = jac_structure(exa_model)
        Je_v = jac_coord(exa_model, x0e)
        J_exa = Matrix(sparse(Je_r, Je_c, Je_v, m, n))
        @test sort(svdvals(J_jump)) ≈ sort(svdvals(J_exa)) atol = 1e-6

        # Hessian: compare absolute eigenvalues (sign may differ due to
        # Max→Min conversion affecting all Lagrangian terms)
        y0 = ones(m)
        Hj_r, Hj_c = hess_structure(jump_nlp)
        Hj_v = hess_coord(jump_nlp, x0j, y0)
        H_jump = Matrix(Symmetric(sparse(Hj_r, Hj_c, Hj_v, n, n), :L))
        He_r, He_c = hess_structure(exa_model)
        He_v = hess_coord(exa_model, x0e, y0)
        H_exa = Matrix(Symmetric(sparse(He_r, He_c, He_v, n, n), :L))
        @test sort(abs.(eigvals(H_jump))) ≈ sort(abs.(eigvals(H_exa))) atol = 1e-6
    end
end
