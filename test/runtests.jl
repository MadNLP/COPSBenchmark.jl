
using Test
using JuMP
using Ipopt
using COPSBenchmark

COPS_INSTANCES = [
    (COPSBenchmark.bearing_model, (50, 50), -1.5482e-1),
    (COPSBenchmark.chain_model, (800,), 5.06891),
    (COPSBenchmark.camshape_model, (1000,), 4.2791), # TODO: result is slightly different
    (COPSBenchmark.catmix_model, (100,), -4.80556e-2),
    (COPSBenchmark.channel_model, (200,), 1.0),
    (COPSBenchmark.elec_model, (50,), 1.0552e3),
    (COPSBenchmark.gasoil_model, (100,), 5.2366e-3),
    (COPSBenchmark.glider_model, (100,), 1.25505e3),
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
    (COPSBenchmark.henon_model, (10,), 6.667736), # N.B: objective depends on the optimizer used.
    (COPSBenchmark.lane_emden_model, (20,), 9.11000),
    (COPSBenchmark.triangle_deer_model, (), 2.01174e3),
    (COPSBenchmark.triangle_pacman_model, (), 1.25045e3),
    (COPSBenchmark.triangle_turtle_model, (), 4.21523e3),
]

@testset "Instance $instance" for (instance, params, result) in COPS_INSTANCES
    model = instance(params...)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    @test JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    # Test that the objective matches the value reported in http://www.mcs.anl.gov/~more/cops/cops3.pdf
    @test JuMP.objective_value(model) ≈ result rtol = 1e-4
end
