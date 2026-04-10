# Electrons on a Sphere Problem.
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.elec_model(::ExaModelsBackend, np; seed = 2713, T = Float64, backend = nothing, kwargs...)
    Random.seed!(seed)

    theta = (2pi) .* rand(np)
    phi = pi .* rand(np)

    core = ExaModels.ExaCore(T; backend = backend)
    ExaModels.@add_variable(core, x, 1:np; start = [cos(theta[i])*sin(phi[i]) for i=1:np])
    ExaModels.@add_variable(core, y, 1:np; start = [sin(theta[i])*sin(phi[i]) for i=1:np])
    ExaModels.@add_variable(core, z, 1:np; start = [cos(phi[i]) for i=1:np])

    itr = [(i,j) for i in 1:np-1 for j in i+1:np]
    ExaModels.@add_objective(core, COPSBenchmark.elec_coulomb_potential(x[i], y[i], z[i], x[j], y[j], z[j]) for (i,j) in itr)

    ExaModels.@add_constraint(core, c1, x[i]^2 + y[i]^2 + z[i]^2 - 1 for i=1:np)

    return ExaModels.ExaModel(core; kwargs...)
end
