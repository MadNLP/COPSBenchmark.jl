# Minimize the sum of the inverse weighted mean ratio of the elements in a fixed-boundary
# tetrahedral mesh by adjusting the locations of the free vertices.

#  This is problem 19 in the COPS (Version 3) collection of
#   E. Dolan and J. More
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

include(joinpath(@__DIR__, "..", "..", "data", "tetra.jl"))

@inline function _tetra_model_exa(x0, TETS::Vector{Int64}, Const::Vector{Int64}; T = Float64, backend = nothing, kwargs...)
    τ = 0.0
    n = length(x0)
    N = round(Int, n / 3)
    E = round(Int, length(TETS) / 4)

    lvar = fill(T(-Inf), n)
    lvar[Const] .= x0[Const]
    lvar[Const .+ N] .= x0[Const .+ N]
    lvar[Const .+ 2 * N] .= x0[Const .+ 2 * N]

    uvar = fill(T(Inf), n)
    uvar[Const] .= x0[Const]
    uvar[Const .+ N] .= x0[Const .+ N]
    uvar[Const .+ 2 * N] .= x0[Const .+ 2 * N]

    c = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(c, x, n; lvar = lvar, uvar = uvar, start = x0)

    # Precompute iteration data: (e, t0, t1, t2, t3) for each element
    itr = [(e, TETS[e], TETS[e + E], TETS[e + 2 * E], TETS[e + 3 * E]) for e in 1:E]

    # Objective: sum over elements of numerator / denominator
    # numerator = sum_i=0:2 ( (x[t1+N*i]-x[t0+N*i])^2 + (2*x[t2+N*i]-x[t1+N*i]-x[t0+N*i])^2/3
    #                        + (3*x[t3+N*i]-x[t2+N*i]-x[t1+N*i]-x[t0+N*i])^2/6 )
    # denominator = 3 * ( sum_i=0:2 (x[t1+N*i]-x[t0+N*i]) * cross_component * sqrt(2) )^(2/3)
    ExaModels.@add_objective(c,
        (
            (x[t1] - x[t0])^2 + (2*x[t2] - x[t1] - x[t0])^2/3 + (3*x[t3] - x[t2] - x[t1] - x[t0])^2/6 +
            (x[t1+N] - x[t0+N])^2 + (2*x[t2+N] - x[t1+N] - x[t0+N])^2/3 + (3*x[t3+N] - x[t2+N] - x[t1+N] - x[t0+N])^2/6 +
            (x[t1+2*N] - x[t0+2*N])^2 + (2*x[t2+2*N] - x[t1+2*N] - x[t0+2*N])^2/3 + (3*x[t3+2*N] - x[t2+2*N] - x[t1+2*N] - x[t0+2*N])^2/6
        ) / (
            3 * (
                (x[t1] - x[t0]) * ((x[t2+N] - x[t0+N])*(x[t3+2*N] - x[t0+2*N]) - (x[t2+2*N] - x[t0+2*N])*(x[t3+N] - x[t0+N])) * sqrt(2) +
                (x[t1+N] - x[t0+N]) * ((x[t2+2*N] - x[t0+2*N])*(x[t3] - x[t0]) - (x[t2] - x[t0])*(x[t3+2*N] - x[t0+2*N])) * sqrt(2) +
                (x[t1+2*N] - x[t0+2*N]) * ((x[t2] - x[t0])*(x[t3+N] - x[t0+N]) - (x[t2+N] - x[t0+N])*(x[t3] - x[t0])) * sqrt(2)
            )^(2/3)
        )
        for (e, t0, t1, t2, t3) in itr
    )

    # Constraint: volume >= τ for each element
    ExaModels.@add_constraint(c, c1,
        (x[t1] - x[t0]) * ((x[t2+N] - x[t0+N])*(x[t3+2*N] - x[t0+2*N]) - (x[t2+2*N] - x[t0+2*N])*(x[t3+N] - x[t0+N])) * sqrt(2) +
        (x[t1+N] - x[t0+N]) * ((x[t2+2*N] - x[t0+2*N])*(x[t3] - x[t0]) - (x[t2] - x[t0])*(x[t3+2*N] - x[t0+2*N])) * sqrt(2) +
        (x[t1+2*N] - x[t0+2*N]) * ((x[t2] - x[t0])*(x[t3+N] - x[t0+N]) - (x[t2+N] - x[t0+N])*(x[t3] - x[t0])) * sqrt(2)
        for (e, t0, t1, t2, t3) in itr;
        lcon = τ, ucon = Inf
    )

    return ExaModels.ExaModel(c; kwargs...)
end

COPSBenchmark.tetra_duct12_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_duct12, TETS_duct12, Const_duct12; kwargs...)
COPSBenchmark.tetra_duct15_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_duct15, TETS_duct15, Const_duct15; kwargs...)
COPSBenchmark.tetra_duct20_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_duct20, TETS_duct20, Const_duct20; kwargs...)
COPSBenchmark.tetra_hook_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_hook, TETS_hook, Const_hook; kwargs...)
COPSBenchmark.tetra_foam5_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_foam5, TETS_foam5, Const_foam5; kwargs...)
COPSBenchmark.tetra_gear_model(::ExaModelsBackend; kwargs...) = _tetra_model_exa(xe_gear, TETS_gear, Const_gear; kwargs...)
