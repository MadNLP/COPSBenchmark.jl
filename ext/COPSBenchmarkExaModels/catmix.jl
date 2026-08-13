# Catalyst Mixing Problem
# Collocation formulation
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# The uniform bounds and zero starts were written as full arrays sized by nh;
# as scalars they say the same thing without depending on the size.  The two
# patterned starts do depend on it and travel as data.
@inline function COPSBenchmark.catmix_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    ne = 2
    nc = 3

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

    core, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.catmix_data, nh)

    h = one(T) / nh   # Final time / nh

    ExaModels.@add_var(core, u, nh, nc; lvar = zero(T), uvar = one(T), start = zero(T))
    ExaModels.@add_var(core, v, nh, ne; start = d.v_start)
    ExaModels.@add_var(core, w, nh, nc, ne; start = zero(T))
    ExaModels.@add_var(core, pp, nh, nc, ne; start = d.pp_start)
    ExaModels.@add_var(core, Dpp, nh, nc, ne; start = zero(T))
    ExaModels.@add_var(core, ppf, ne; start = T[mod(i,ne) for i in 1:ne])

    ExaModels.@add_obj(core, neg_one + ppf[1] + ppf[2])
    ExaModels.@add_obj(core, alpha/h*(u[i+1, j] - u[i, j])^2 for i in 1:nh-1, j in 1:nc)

    ExaModels.@add_con(
        core,
        c1,
        pp[i, k, s] - v[i, s] - h*sum(w[i, j, s]*(rho^j/T(factorial(j))) for j in 1:nc) for i=1:nh,  (k, rho) in rho_index, s=1:ne
    )

    ExaModels.@add_con(
        core,
        c2,
        Dpp[i, k, s] - sum(w[i, j, s]*(rho^(j-1)/T(factorial(j-1))) for j in 1:nc) for i=1:nh, (k, rho) in rho_index, s=1:ne
    )

    ExaModels.@add_con(
        core,
        c3,
        ppf[s] - v[k, s] - h * sum(w[k, j, s] / T(factorial(j)) for j in 1:nc) for k in nh:nh, s in 1:ne
    )

    ExaModels.@add_con(
        core,
        c4,
        v[i, s] + sum(w[i, j, s] * h / T(factorial(j)) for j in 1:nc) - v[i+1, s] for i in 1:nh-1, s in 1:ne
    )

    ExaModels.@add_con(
        core,
        c5,
        Dpp[i,j,1] - u[i,j] * (ten_T*pp[i,j,2] - pp[i,j,1]) for i=1:nh, j=1:nc
    )

    ExaModels.@add_con(
        core,
        c6,
        Dpp[i,j,2] - u[i,j] * (pp[i,j,1] - ten_T*pp[i,j,2]) + (one_T - u[i,j])*pp[i,j,2] for i=1:nh, j=1:nc
    )

    ExaModels.@add_con(
        core,
        c7,
        v[1, s] - bc for (s, bc) in [(i, bc[i]) for i in 1:ne]
    )

    return core
end

@inline COPSBenchmark.catmix_args(::ExaModelsBackend, nh; T = Float64) = (nh,)

@inline COPSBenchmark.catmix_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.catmix_recipe(b; T = T, backend = backend),
        COPSBenchmark.catmix_args(b, nh; T = T)...;
        kwargs...,
    )
