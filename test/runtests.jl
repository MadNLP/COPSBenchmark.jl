
using Test
using JuMP
using Ipopt
using ExaModels
using NLPModels
using NLPModelsIpopt
using NLPModelsJuMP
using SparseArrays
using LinearAlgebra
using COPSBenchmark

COPS_INSTANCES = [
    (COPSBenchmark.bearing_model, (50, 50), -1.5482e-1),
    (COPSBenchmark.chain_model, (800,), 5.06891),
    (COPSBenchmark.camshape_model, (1000,), 4.2791), # TODO: result is slightly different
    (COPSBenchmark.catmix_model, (100,), -4.80556e-2),
    (COPSBenchmark.channel_model, (200,), 1.0),
    (COPSBenchmark.elec_model, (50,), 1.0552e3),
    (COPSBenchmark.gasoil_model, (100,), 5.2366e-3),
    (COPSBenchmark.glider_model, (100,), -1.25505e3),
    (COPSBenchmark.marine_model, (100,), 1.97462e7),
    (COPSBenchmark.methanol_model, (100,), 9.02229e-3),
    (COPSBenchmark.minsurf_model, (50, 50), 2.51488),
    (COPSBenchmark.minsurf_model, (50, 75), 2.50568),
    (COPSBenchmark.minsurf_model, (50, 100), 2.50694),
    (COPSBenchmark.pinene_model, (100,), 1.98721e1),
    (COPSBenchmark.polygon_model, (100,), -0.674981), # N.B: objective depends on the optimizer used.
    (COPSBenchmark.robot_model, (200,), 9.14138),
    (COPSBenchmark.rocket_model, (400,), 1.01283),
    (COPSBenchmark.steering_model, (200,), 5.54577e-1),
    (COPSBenchmark.tetra_duct15_model, (), 1.04951e4),
    (COPSBenchmark.tetra_duct20_model, (), 4.82685e3),
    (COPSBenchmark.tetra_foam5_model, (), 6.42560e3),
    (COPSBenchmark.tetra_gear_model, (), 4.15163e3),
    (COPSBenchmark.tetra_hook_model, (), 6.05735e3),
    (COPSBenchmark.torsion_model, (50, 50), -4.18087e-1),
    (COPSBenchmark.dirichlet_model, (20,), 1.71464e-2),
    # N.B.: comment henon and lane-emden as they take too long in the CI.
    # (COPSBenchmark.henon_model, (10,), 6.667736), # N.B: objective depends on the optimizer used.
    # (COPSBenchmark.lane_emden_model, (20,), 9.11000),
    (COPSBenchmark.triangle_deer_model, (), 2.01174e3),
    (COPSBenchmark.triangle_pacman_model, (), 1.25045e3),
    (COPSBenchmark.triangle_turtle_model, (), 4.21523e3),
]

function solve_backend(model, ::COPSBenchmark.JuMPBackend)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    obj_val = JuMP.objective_value(model)
    return status, obj_val
end

function solve_backend(model, ::COPSBenchmark.ExaModelsBackend)
    results = ipopt(model; print_level=0, tol=1e-8)
    status = results.status == :first_order
    obj_val = results.objective
    return status, obj_val
end

@testset "Backend $(backend)" for backend in [
    COPSBenchmark.JuMPBackend(),
    COPSBenchmark.ExaModelsBackend(),
]
    @testset "Instance $instance" for (instance, params, result) in COPS_INSTANCES
        if hasmethod(instance, Tuple{typeof.(params)..., typeof(backend)})
            model = instance(params..., backend)
            status, obj_val = solve_backend(model, backend)
            # Test that the objective matches the value reported in
            # http://www.mcs.anl.gov/~more/cops/cops3.pdf
            @test obj_val ≈ result rtol = 1e-4
            @test status
        end
    end
end

# Models without lifted variables (same variable/constraint space in both backends).
COMPARE_INSTANCES = [
    (COPSBenchmark.bearing_model, (10, 10)),
    (COPSBenchmark.camshape_model, (100,)),
    (COPSBenchmark.catmix_model, (10,)),
    (COPSBenchmark.chain_model, (100,)),
    (COPSBenchmark.elec_model, (25,)),
    (COPSBenchmark.gasoil_model, (10,)),
    (COPSBenchmark.marine_model, (10,)),
    (COPSBenchmark.methanol_model, (10,)),
    (COPSBenchmark.minsurf_model, (10, 10)),
    (COPSBenchmark.pinene_model, (10,)),
    (COPSBenchmark.rocket_model, (50,)),
    (COPSBenchmark.steering_model, (50,)),
]

@testset "Compare callbacks: $instance" for (instance, params) in COMPARE_INSTANCES
    jump_model = instance(params..., COPSBenchmark.JuMPBackend())
    exa_model = instance(params..., COPSBenchmark.ExaModelsBackend())

    jump_nlp = MathOptNLPModel(jump_model)

    # Structural comparison
    @test get_nvar(jump_nlp) == get_nvar(exa_model)
    @test get_ncon(jump_nlp) == get_ncon(exa_model)

    # Variable bounds
    @test get_lvar(jump_nlp) ≈ get_lvar(exa_model)
    @test get_uvar(jump_nlp) ≈ get_uvar(exa_model)

    # Constraint bounds
    @test get_lcon(jump_nlp) ≈ get_lcon(exa_model)
    @test get_ucon(jump_nlp) ≈ get_ucon(exa_model)

    # Starting point
    @test jump_nlp.meta.x0 ≈ exa_model.meta.x0

    # Callbacks at starting point
    x0 = jump_nlp.meta.x0
    n = get_nvar(jump_nlp)
    m = get_ncon(jump_nlp)
    @test obj(jump_nlp, x0) ≈ obj(exa_model, x0) rtol = 1e-6
    @test grad(jump_nlp, x0) ≈ grad(exa_model, x0) rtol = 1e-6
    @test cons(jump_nlp, x0) ≈ cons(exa_model, x0) rtol = 1e-6

    # Jacobian structure and values (ordering may differ between backends)
    Jj_r, Jj_c = jac_structure(jump_nlp)
    Jj_v = jac_coord(jump_nlp, x0)
    J_jump = sparse(Jj_r, Jj_c, Jj_v, m, n)

    Je_r, Je_c = jac_structure(exa_model)
    Je_v = jac_coord(exa_model, x0)
    J_exa = sparse(Je_r, Je_c, Je_v, m, n)

    @test sort(collect(zip(Jj_r, Jj_c))) == sort(collect(zip(Je_r, Je_c)))
    @test Matrix(J_jump) ≈ Matrix(J_exa) rtol = 1e-6

    # Hessian structure and values (ordering may differ between backends)
    y0 = ones(m)
    Hj_r, Hj_c = hess_structure(jump_nlp)
    Hj_v = hess_coord(jump_nlp, x0, y0)
    H_jump = Symmetric(sparse(Hj_r, Hj_c, Hj_v, n, n), :L)

    He_r, He_c = hess_structure(exa_model)
    He_v = hess_coord(exa_model, x0, y0)
    H_exa = Symmetric(sparse(He_r, He_c, He_v, n, n), :L)

    @test sort(collect(zip(Hj_r, Hj_c))) == sort(collect(zip(He_r, He_c)))
    @test Matrix(H_jump) ≈ Matrix(H_exa) rtol = 1e-6
end
