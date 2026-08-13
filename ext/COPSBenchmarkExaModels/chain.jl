# Hanging Chain

# Find the chain (of uniform density) of length L suspended between two points with minimal
# potential energy.

#   This is problem 4 in the COPS (Version 3) collection of
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

# The public size `n` is not the discretization: `nh = max(2, div(n-4, 4))` is.
# `chain_args` does that reduction and builds the four starting curves, which
# are comprehensions over nh; the recipe is written against nh directly.
@inline function COPSBenchmark.chain_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    L = 4
    a = 1
    b = 3
    tf = 1.0

    c, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.chain_data, nh)

    h = tf / nh

    ExaModels.@add_var(c, u, nh + 1; start = d.u0)
    ExaModels.@add_var(c, x1, nh + 1; start = d.x10)
    ExaModels.@add_var(c, x2, nh + 1; start = d.x20)
    ExaModels.@add_var(c, x3, nh + 1; start = d.x30)

    # Indexes the last point, so the index arrives as data.
    ExaModels.@add_obj(c, x2[k] for k in (nh+1):(nh+1))

    ExaModels.@add_con(
        c,
        c1,
        x1[j + 1] - x1[j] - 1 / 2 * h * (u[j] + u[j + 1]) for j in 1:nh
    )

    ExaModels.@add_con(
        c,
        c2,
        x1[1] - a
    )

    ExaModels.@add_con(
        c,
        c3,
        x1[k] - b for k in (nh+1):(nh+1)
    )

    ExaModels.@add_con(
        c,
        c4,
        x2[1]
    )

    ExaModels.@add_con(
        c,
        c5,
        x3[1]
    )

    ExaModels.@add_con(
        c,
        c6,
        x3[k] - L for k in (nh+1):(nh+1)
    )

    ExaModels.@add_con(
        c,
        c7,
        x2[j + 1] - x2[j] - 1 / 2 * h * (x1[j] * sqrt(1 + u[j]^2) + x1[j + 1] * sqrt(1 + u[j + 1]^2)) for j in 1:nh
    )

    ExaModels.@add_con(
        c,
        c8,
        x3[j + 1] - x3[j] - 1 / 2 * h * (sqrt(1 + u[j]^2) + sqrt(1 + u[j + 1]^2)) for j in 1:nh
    )

    return c
end

@inline COPSBenchmark.chain_args(::ExaModelsBackend, n; T = Float64) = (max(2, div(n - 4, 4)),)

@inline COPSBenchmark.chain_model(b::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.chain_recipe(b; T = T, backend = backend),
        COPSBenchmark.chain_args(b, n; T = T)...;
        kwargs...,
    )
