# Find the polygon of maximal area, among polygons with nv sides and
#   diameter d <= 1

#   This is problem 1 in the COPS (Version 3) collection of
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

# The public size `n` is halved to get the vertex count, which `polygon_args`
# does along with the angle start and the diameter pair list -- both
# comprehensions over that count.  The two rows fixing the last vertex index it,
# so they come from one-element sets.
@inline function COPSBenchmark.polygon_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    c, N = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.polygon_data, N)

    ExaModels.@add_var(c, r, N; lvar = 0.0, uvar = 1.0, start = 1.0)
    ExaModels.@add_var(c, θ, N; lvar = 0.0, uvar = T(π), start = d.θ0)

    # Objective: maximize area = 0.5 * sum(r[i]*r[i+1]*sin(θ[i+1]-θ[i]))
    ExaModels.@add_obj(c, -0.5 * r[i] * r[i+1] * sin(θ[i+1] - θ[i]) for i in 1:N-1)

    # Fix last angle and radius
    ExaModels.@add_con(c, c1, θ[k] - T(π) for k in N:N)
    ExaModels.@add_con(c, c2, r[k] for k in N:N)

    # Impose ordering on angles: θ[i+1] >= θ[i]
    ExaModels.@add_con(c, c3, θ[i+1] - θ[i] for i in 1:N-1; lcon = 0.0, ucon = Inf)

    # Diameter constraint: r[i]^2 + r[j]^2 - 2*r[i]*r[j]*cos(θ[i]-θ[j]) <= 1
    ExaModels.@add_con(c, c4, r[i]^2 + r[j]^2 - 2*r[i]*r[j]*cos(θ[i] - θ[j]) - 1 for (i, j) in d.pairs; lcon = -Inf, ucon = 0.0)

    return c
end

@inline COPSBenchmark.polygon_args(::ExaModelsBackend, n::Int) = (div(n, 2),)

@inline COPSBenchmark.polygon_model(b::ExaModelsBackend, n::Int; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.polygon_recipe(b; T = T, backend = backend),
        COPSBenchmark.polygon_args(b, n)...;
        kwargs...,
    )
