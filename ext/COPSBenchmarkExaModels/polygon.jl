# Find the polygon of maximal area, among polygons with nv sides and
#   diameter d <= 1

#   This is problem 1 in the COPS (Version 3) collection of
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

@inline function COPSBenchmark.polygon_model(::ExaModelsBackend, n::Int; T = Float64, backend = nothing, kwargs...)
    N = div(n, 2)

    c = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@var(c, r, N; lvar = 0.0, uvar = 1.0, start = 1.0)
    ExaModels.@var(c, θ, N; lvar = 0.0, uvar = T(π), start = [i * π / (N - 1) - π / (N - 1) for i in 1:N])

    # Objective: maximize area = 0.5 * sum(r[i]*r[i+1]*sin(θ[i+1]-θ[i]))
    ExaModels.@obj(c, -0.5 * r[i] * r[i+1] * sin(θ[i+1] - θ[i]) for i in 1:N-1)

    # Fix last angle and radius
    ExaModels.@con(c, c1, θ[N] - T(π))
    ExaModels.@con(c, c2, r[N])

    # Impose ordering on angles: θ[i+1] >= θ[i]
    ExaModels.@con(c, c3, θ[i+1] - θ[i] for i in 1:N-1; lcon = 0.0, ucon = Inf)

    # Diameter constraint: r[i]^2 + r[j]^2 - 2*r[i]*r[j]*cos(θ[i]-θ[j]) <= 1
    pairs = [(i, j) for i in 1:N-1 for j in i+1:N]
    ExaModels.@con(c, c4, r[i]^2 + r[j]^2 - 2*r[i]*r[j]*cos(θ[i] - θ[j]) - 1 for (i, j) in pairs; lcon = -Inf, ucon = 0.0)

    return ExaModels.ExaModel(c; kwargs...)
end
