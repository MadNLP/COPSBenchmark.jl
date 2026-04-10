# Hanging Chain
# Find the chain (of uniform density) of length L suspended between two points
# with minimal potential energy.

@inline function COPSBenchmark.chain_model(::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...)
    nh    = COPSBenchmark.chain_nh(n)
    L     = COPSBenchmark.chain_L
    a     = COPSBenchmark.chain_a
    b     = COPSBenchmark.chain_b
    h     = COPSBenchmark.chain_tf / nh

    c = ExaModels.ExaCore(T; backend = backend)
    ExaModels.@add_variable(c, u, nh + 1; start = [COPSBenchmark.chain_u_start(k, nh) for k in 1:nh+1])
    ExaModels.@add_variable(c, x1, nh + 1; start = [COPSBenchmark.chain_x1_start(k, nh) for k in 1:nh+1])
    ExaModels.@add_variable(c, x2, nh + 1; start = [COPSBenchmark.chain_x2_start(k, nh) for k in 1:nh+1])
    ExaModels.@add_variable(c, x3, nh + 1; start = [COPSBenchmark.chain_u_start(k, nh) for k in 1:nh+1])

    ExaModels.@add_objective(c, x2[nh + 1])

    ExaModels.@add_constraint(c, c1, x1[j + 1] - x1[j] - 1/2 * h * (u[j] + u[j + 1]) for j in 1:nh)
    ExaModels.@add_constraint(c, c2, x1[1] - a)
    ExaModels.@add_constraint(c, c3, x1[nh + 1] - b)
    ExaModels.@add_constraint(c, c4, x2[1])
    ExaModels.@add_constraint(c, c5, x3[1])
    ExaModels.@add_constraint(c, c6, x3[nh+1] - L)
    ExaModels.@add_constraint(c, c7, x2[j + 1] - x2[j] - 1/2 * h * (x1[j] * sqrt(1 + u[j]^2) + x1[j + 1] * sqrt(1 + u[j + 1]^2)) for j in 1:nh)
    ExaModels.@add_constraint(c, c8, x3[j + 1] - x3[j] - 1/2 * h * (sqrt(1 + u[j]^2) + sqrt(1 + u[j + 1]^2)) for j in 1:nh)

    return ExaModels.ExaModel(c; kwargs...)
end
