# Electrons on a Sphere Problem.
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.elec_model(np, ::ExaModelsBackend; seed = 2713, T = Float64, backend = nothing, kwargs...)
    Random.seed!(seed)

    # Set the starting point to a quasi-uniform distribution
    # of electrons on a unit sphere
    theta = (2pi) .* rand(np)
    phi = pi .* rand(np)

    core = ExaModels.ExaCore(T; backend= backend)
    x = ExaModels.variable(core, 1:np; start = [cos(theta[i])*sin(phi[i]) for i=1:np])
    y = ExaModels.variable(core, 1:np; start = [sin(theta[i])*sin(phi[i]) for i=1:np])
    z = ExaModels.variable(core, 1:np; start = [cos(phi[i]) for i=1:np])

    # Coulomb potential
    itr = [(i,j) for i in 1:np-1 for j in i+1:np]
    ExaModels.objective(core, 1.0 / sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2) for (i,j) in itr)

    # Unit-ball
    ExaModels.constraint(core, x[i]^2 + y[i]^2 + z[i]^2 - 1 for i=1:np)

    return ExaModels.ExaModel(core; kwargs...)
end
