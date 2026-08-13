# Electrons on a Sphere Problem.
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# The quasi-uniform starting distribution is a random draw, and the Coulomb
# pair list is a double comprehension over the size.  Both are data: they are
# drawn and built in `elec_args`, which owns the seed, so the recipe holds no
# randomness and the structure is the same whatever was drawn.
@inline function COPSBenchmark.elec_recipe(
    ::ExaModelsBackend; seed = 2713, T = Float64, backend = nothing,
)
    core, np = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgCall(COPSBenchmark.elec_data, (np, seed))

    ExaModels.@add_var(core, x, 1:np; start = d.x0)
    ExaModels.@add_var(core, y, 1:np; start = d.y0)
    ExaModels.@add_var(core, z, 1:np; start = d.z0)

    # Coulomb potential
    ExaModels.@add_obj(core, 1.0 / sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2) for (i,j) in d.itr)

    # Unit-ball
    ExaModels.@add_con(core, c1, x[i]^2 + y[i]^2 + z[i]^2 - 1 for i=1:np)

    return core
end

@inline COPSBenchmark.elec_args(::ExaModelsBackend, np) = (np,)

@inline COPSBenchmark.elec_model(b::ExaModelsBackend, np; seed = 2713, T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.elec_recipe(b; seed = seed, T = T, backend = backend),
        COPSBenchmark.elec_args(b, np)...;
        kwargs...,
    )
