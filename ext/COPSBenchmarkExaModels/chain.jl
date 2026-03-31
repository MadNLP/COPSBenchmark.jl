# Hanging Chain

# Find the chain (of uniform density) of length L suspended between two points with minimal
# potential energy.

#   This is problem 4 in the COPS (Version 3) collection of
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

function COPSBenchmark.chain_model(::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...)
    nh = max(2, div(n - 4, 4))

    L = 4
    a = 1
    b = 3
    tmin = b > a ? 1 / 4 : 3 / 4
    tf = 1.0
    h = tf / nh

    c = ExaModels.ExaCore(T; backend = backend)
    ExaModels.@var(c, u, nh + 1; start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])
    ExaModels.@var(c, x1, nh + 1; start = [4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a for k in 1:nh+1])
    ExaModels.@var(c, x2, nh + 1; start = [(4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a) *
        (4 * abs(b - a) * (k / nh - tmin)) for k in 1:nh+1])
    ExaModels.@var(c, x3, nh + 1;  start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])

    ExaModels.@obj(c, x2[nh + 1])

    ExaModels.@con(
        c,
        c1,
        x1[j + 1] - x1[j] - 1 / 2 * h * (u[j] + u[j + 1]) for j in 1:nh
    )

    ExaModels.@con(
        c,
        c2,
        x1[1] - a
    )

    ExaModels.@con(
        c,
        c3,
        x1[nh + 1] - b
    )

    ExaModels.@con(
        c,
        c4,
        x2[1]
    )

    ExaModels.@con(
        c,
        c5,
        x3[1]
    )

    ExaModels.@con(
        c,
        c6,
        x3[nh+1] - L
    )

    ExaModels.@con(
        c,
        c7,
        x2[j + 1] - x2[j] - 1 / 2 * h * (x1[j] * sqrt(1 + u[j]^2) + x1[j + 1] * sqrt(1 + u[j + 1]^2)) for j in 1:nh
    )

    ExaModels.@con(
        c,
        c8,
        x3[j + 1] - x3[j] - 1 / 2 * h * (sqrt(1 + u[j]^2) + sqrt(1 + u[j + 1]^2)) for j in 1:nh
    )

    return ExaModels.ExaModel(c; kwargs...)
end

