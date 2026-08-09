# Minimize the sum of the inverse weighted mean ratio of the elements in a fixed-boundary
# triangular mesh by adjusting the locations of the free vertices.

#  This is problem 18 in the COPS (Version 3) collection of
#   E. Dolan and J. More
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

include(joinpath(@__DIR__, "..", "..", "data", "triangle.jl"))

@inline function _triangle_model_exa(x0, TRIS::Vector{Int64}, Const::Vector{Int64}; T = Float64, backend = nothing, kwargs...)
    τ = 0.0
    n = length(x0)
    N = Int(div(n, 2))
    E = Int(div(length(TRIS), 3))

    lvar = fill(T(-Inf), n)
    lvar[Const] = x0[Const]
    lvar[Const .+ N] = x0[Const .+ N]

    uvar = fill(T(Inf), n)
    uvar[Const] = x0[Const]
    uvar[Const .+ N] = x0[Const .+ N]

    c = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))

    ExaModels.@add_var(c, x, n; lvar = lvar, uvar = uvar, start = x0)

    # Precompute iteration data: (e, t0, t1, t2) for each element
    itr = [(e, TRIS[e], TRIS[e + E], TRIS[e + 2 * E]) for e in 1:E]

    # Objective: sum over elements of (numerator / denominator)
    # numerator = sum_i ((x[t1+N*i] - x[t0+N*i])^2 + (2*x[t2+N*i] - x[t1+N*i] - x[t0+N*i])^2/3 for i=0:1)
    # denominator = 2 * (2 * ((x[t1]-x[t0])*(x[t2+N]-x[t0+N]) - (x[t2]-x[t0])*(x[t1+N]-x[t0+N])) / sqrt(3))
    ExaModels.@add_obj(c,
        ((x[t1] - x[t0])^2 + (2*x[t2] - x[t1] - x[t0])^2/3 +
         (x[t1+N] - x[t0+N])^2 + (2*x[t2+N] - x[t1+N] - x[t0+N])^2/3) /
        (2 * 2 * ((x[t1]-x[t0])*(x[t2+N]-x[t0+N]) - (x[t2]-x[t0])*(x[t1+N]-x[t0+N])) / sqrt(3))
        for (e, t0, t1, t2) in itr
    )

    # Constraint: area >= τ for each element
    ExaModels.@add_con(c, c1,
        2 * ((x[t1]-x[t0])*(x[t2+N]-x[t0+N]) - (x[t2]-x[t0])*(x[t1+N]-x[t0+N])) / sqrt(3)
        for (e, t0, t1, t2) in itr;
        lcon = τ, ucon = Inf
    )

    return ExaModels.ExaModel(c; kwargs...)
end

COPSBenchmark.triangle_deer_model(::ExaModelsBackend; kwargs...) = _triangle_model_exa(xe_deer, TRIS_deer, Const_deer; kwargs...)
COPSBenchmark.triangle_pacman_model(::ExaModelsBackend; kwargs...) = _triangle_model_exa(xe_pacman, TRIS_pacman, Const_pacman; kwargs...)
COPSBenchmark.triangle_turtle_model(::ExaModelsBackend; kwargs...) = _triangle_model_exa(xe_turtle, TRIS_turtle, Const_turtle; kwargs...)

# 1:1 transcription of _triangle_model_exa (triangle.jl) returning the core.
function _lazy_triangle(x0, TRIS, Const)
    τ = 0.0
    n = length(x0)
    N = Int(div(n, 2))
    E = Int(div(length(TRIS), 3))

    lvar = fill(-Inf, n)
    lvar[Const] = x0[Const]
    lvar[Const .+ N] = x0[Const .+ N]

    uvar = fill(Inf, n)
    uvar[Const] = x0[Const]
    uvar[Const .+ N] = x0[Const .+ N]

    c = ExaCore(concrete = Val(true))
    c, x = add_var(c, n; lvar = lvar, uvar = uvar, start = x0)

    itr = [(e, TRIS[e], TRIS[e + E], TRIS[e + 2 * E]) for e in 1:E]

    c, _ = add_obj(c,
        ((x[t1] - x[t0])^2 + (2*x[t2] - x[t1] - x[t0])^2/3 +
         (x[t1+N] - x[t0+N])^2 + (2*x[t2+N] - x[t1+N] - x[t0+N])^2/3) /
        (2 * 2 * ((x[t1]-x[t0])*(x[t2+N]-x[t0+N]) - (x[t2]-x[t0])*(x[t1+N]-x[t0+N])) / sqrt(3))
        for (e, t0, t1, t2) in itr)

    c, _ = add_con(c,
        (2 * ((x[t1]-x[t0])*(x[t2+N]-x[t0+N]) - (x[t2]-x[t0])*(x[t1+N]-x[t0+N])) / sqrt(3)
         for (e, t0, t1, t2) in itr);
        lcon = τ, ucon = Inf)
    c
end

COPSBenchmark.triangle_deer_core() = _lazy_triangle(xe_deer, TRIS_deer, Const_deer)
COPSBenchmark.triangle_pacman_core() = _lazy_triangle(xe_pacman, TRIS_pacman, Const_pacman)
COPSBenchmark.triangle_turtle_core() = _lazy_triangle(xe_turtle, TRIS_turtle, Const_turtle)
