# Hanging Chain

# Find the chain (of uniform density) of length L suspended between two points with minimal
# potential energy.

#   This is problem 4 in the COPS (Version 3) collection of
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

@inline function COPSBenchmark.chain_model(::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...)
    nh = max(2, div(n - 4, 4))

    L = 4
    a = 1
    b = 3
    tmin = b > a ? 1 / 4 : 3 / 4
    tf = 1.0
    h = tf / nh

    c = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))
    ExaModels.@add_var(c, u, nh + 1; start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])
    ExaModels.@add_var(c, x1, nh + 1; start = [4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a for k in 1:nh+1])
    ExaModels.@add_var(c, x2, nh + 1; start = [(4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a) *
        (4 * abs(b - a) * (k / nh - tmin)) for k in 1:nh+1])
    ExaModels.@add_var(c, x3, nh + 1;  start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])

    ExaModels.@add_obj(c, x2[nh + 1])

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
        x1[nh + 1] - b
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
        x3[nh+1] - L
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

    return ExaModels.ExaModel(c; kwargs...)
end


function COPSBenchmark.chain_core()
    args = ExaModels.ArgTracer()
    c = ExaCore(concrete = Val(true))
    nh = max(2, div(args.n - 4, 4))

    L = 4
    a = 1
    b = 3
    tmin = b > a ? 1 / 4 : 3 / 4
    tf = 1.0
    h = tf / nh

    c, u = add_var(c, nh + 1; start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])
    c, x1 = add_var(c, nh + 1; start = [4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a for k in 1:nh+1])
    c, x2 = add_var(c, nh + 1; start = [(4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a) *
        (4 * abs(b - a) * (k / nh - tmin)) for k in 1:nh+1])
    c, x3 = add_var(c, nh + 1; start = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1])

    c, _ = add_obj(c, x2[m] for m in nh+1:nh+1)

    c, _ = add_con(c, x1[j+1] - x1[j] - 1 / 2 * h * (u[j] + u[j+1]) for j in 1:nh)
    c, _ = add_con(c, x1[1] - a for _ in 1:1)
    c, _ = add_con(c, x1[m] - b for m in nh+1:nh+1)
    c, _ = add_con(c, x2[1] for _ in 1:1)
    c, _ = add_con(c, x3[1] for _ in 1:1)
    c, _ = add_con(c, x3[m] - L for m in nh+1:nh+1)
    c, _ = add_con(c, x2[j+1] - x2[j] - 1 / 2 * h * (x1[j] * sqrt(1 + u[j]^2) + x1[j+1] * sqrt(1 + u[j+1]^2)) for j in 1:nh)
    c, _ = add_con(c, x3[j+1] - x3[j] - 1 / 2 * h * (sqrt(1 + u[j]^2) + sqrt(1 + u[j+1]^2)) for j in 1:nh)
    c
end
