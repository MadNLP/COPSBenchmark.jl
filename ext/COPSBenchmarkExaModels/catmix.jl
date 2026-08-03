# Catalyst Mixing Problem
# Collocation formulation
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.catmix_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    ne = 2
    nc = 3

    h = T(1) / T(nh)   # Final time / nh (was Int/Int → Float64)

    rho = T[
        0.11270166537926,
        0.50000000000000,
        0.88729833462074,
    ]
    bc = T[1, 0]   # Boundary conditions for x
    alpha = T(0)   # Smoothing parameter
    neg_one = -T(1)
    ten_T = T(10)
    one_T = T(1)
    rho_index = [(i, rho[i]) for i in 1:nc]

    c = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))
    ExaModels.@add_var(c, u, nh, nc; lvar = zeros(T, nh, nc), uvar = ones(T, nh, nc), start = zeros(T, nh, nc))
    ExaModels.@add_var(c, v, nh, ne; start = T[mod(j, ne) for i in 1:nh, j in 1:ne])
    ExaModels.@add_var(c, w, nh, nc, ne; start = zeros(T, nh, nc, ne))
    ExaModels.@add_var(c, pp, nh, nc, ne; start = T[mod(k, ne) for i in 1:nh, j in 1:nc, k in 1:ne])
    ExaModels.@add_var(c, Dpp, nh, nc, ne; start = zeros(T, nh, nc, ne))
    ExaModels.@add_var(c, ppf, ne; start = T[mod(i,ne) for i in 1:ne])

    ExaModels.@add_obj(c, neg_one + ppf[1] + ppf[2])
    ExaModels.@add_obj(c, alpha/h*(u[i+1, j] - u[i, j])^2 for i in 1:nh-1, j in 1:nc)

    ExaModels.@add_con(
        c,
        c1,
        pp[i, k, s] - v[i, s] - h*sum(w[i, j, s]*(rho^j/T(factorial(j))) for j in 1:nc) for i=1:nh,  (k, rho) in rho_index, s=1:ne
    )

    ExaModels.@add_con(
        c,
        c2,
        Dpp[i, k, s] - sum(w[i, j, s]*(rho^(j-1)/T(factorial(j-1))) for j in 1:nc) for i=1:nh, (k, rho) in rho_index, s=1:ne
    )

    ExaModels.@add_con(
        c,
        c3,
        ppf[s] - v[nh, s] - h * sum(w[nh, j, s] / T(factorial(j)) for j in 1:nc) for s in 1:ne
    )

    ExaModels.@add_con(
        c,
        c4,
        v[i, s] + sum(w[i, j, s] * h / T(factorial(j)) for j in 1:nc) - v[i+1, s] for i in 1:nh-1, s in 1:ne
    )



    ExaModels.@add_con(
        c,
        c5,
        Dpp[i,j,1] - u[i,j] * (ten_T*pp[i,j,2] - pp[i,j,1]) for i=1:nh, j=1:nc
    )

    ExaModels.@add_con(
        c,
        c6,
        Dpp[i,j,2] - u[i,j] * (pp[i,j,1] - ten_T*pp[i,j,2]) + (one_T - u[i,j])*pp[i,j,2] for i=1:nh, j=1:nc
    )


    ExaModels.@add_con(
        c,
        c7,
        v[1, s] - bc for (s, bc) in [(i, bc[i]) for i in 1:ne]
    )

    return ExaModels.ExaModel(c; kwargs...)
end
