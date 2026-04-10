# Hanging Chain
# Find the chain (of uniform density) of length L suspended between two points
# with minimal potential energy.

# This file has been adapted from https://github.com/JuliaSmoothOptimizers/OptimizationProblems.jl

function COPSBenchmark.chain_model(::JuMPBackend, n::Int)
    nh    = COPSBenchmark.chain_nh(n)
    L     = COPSBenchmark.chain_L
    a     = COPSBenchmark.chain_a
    b     = COPSBenchmark.chain_b
    h     = COPSBenchmark.chain_tf / nh

    nlp = Model()

    @variable(nlp, u[k = 1:(nh + 1)], start = COPSBenchmark.chain_u_start(k, nh))
    @variable(nlp, x1[k = 1:(nh + 1)], start = COPSBenchmark.chain_x1_start(k, nh))
    @variable(nlp, x2[k = 1:(nh + 1)], start = COPSBenchmark.chain_x2_start(k, nh))
    @variable(nlp, x3[k = 1:(nh + 1)], start = COPSBenchmark.chain_u_start(k, nh))

    @objective(nlp, Min, x2[nh + 1])

    for j = 1:nh
        @constraint(nlp, x1[j + 1] - x1[j] - 1/2 * h * (u[j] + u[j + 1]) == 0)
    end
    @constraint(nlp, x1[1] == a)
    @constraint(nlp, x1[nh + 1] == b)
    @constraint(nlp, x2[1] == 0)
    @constraint(nlp, x3[1] == 0)
    @constraint(nlp, x3[nh + 1] == L)

    @constraint(
        nlp,
        [j = 1:nh],
        x2[j + 1] - x2[j] -
        1/2 * h * (x1[j] * sqrt(1 + u[j]^2) + x1[j + 1] * sqrt(1 + u[j + 1]^2)) == 0
    )
    @constraint(
        nlp,
        [j = 1:nh],
        x3[j + 1] - x3[j] - 1/2 * h * (sqrt(1 + u[j]^2) + sqrt(1 + u[j + 1]^2)) == 0
    )

    return nlp
end
