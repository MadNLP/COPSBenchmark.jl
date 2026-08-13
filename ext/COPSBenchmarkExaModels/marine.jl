# Marine Population Dynamics Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# Same shape as gasoil: `h` stays symbolic, the measurement table and the
# starting values built from it are data.  `v` and `uc` share one start array in
# the original (nc = 1, so the shapes coincide) and share one field here.
@inline function COPSBenchmark.marine_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    nc = 1     # number of collocation points
    ne = 8     # number of differential equations
    nm = 21    # number of measurements

    rho = T[0.5]
    tf = T(10)                                    # tau[nm]
    zero_T = T(0)

    core, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.marine_data, nh)

    h = tf / nh                                   # uniform interval length

    # Growth rates
    ExaModels.@add_var(core, g, 1:ne-1; lvar = zero_T)
    # Mortality rates
    ExaModels.@add_var(core, m, 1:ne; lvar = zero_T)
    ExaModels.@add_var(core, v, 1:nh, 1:ne; start=d.v_start)
    ExaModels.@add_var(core, w, 1:nh, 1:nc, 1:ne; start=zero_T)
    ExaModels.@add_var(core, uc, 1:nh, 1:nc, 1:ne; start=d.v_start)
    ExaModels.@add_var(core, Duc, 1:nh, 1:nc, 1:ne; start=zero_T)

    itr2 = [(j,rho[j]) for j=1:nc]

    ExaModels.@add_obj(core, (v[itauj,s] + sum(w[itauj,k,s]*(tauj-tj)^k/(T(factorial(k))*h^(k-1)) for k in 1:nc) - zjs)^2 for (j,s,itauj,tauj, tj, zjs) in d.itr)

    # Collocation model
    ExaModels.@add_con(
        core,
        c1,
        - uc[i, j, s] + v[i,s] + h*sum(w[i,k,s]*(rhoj^k/T(factorial(k))) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    ExaModels.@add_con(
        core,
        c2,
        - Duc[i, j, s] + sum(w[i,k,s]*(rhoj^(k-1)/T(factorial(k-1))) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    # Continuity
    ExaModels.@add_con(
        core,
        c3,
        v[i, s] + sum(w[i, j, s]*h/T(factorial(j)) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    # Boundary conditions
    ExaModels.@add_con(
        core,
        c4,
        - Duc[i, j, 1] -(m[1]+g[1])*uc[i, j, 1]
        for i=1:nh, j=1:nc
    )
    ExaModels.@add_con(
        core,
        c5,
        - Duc[i, j, ne] + g[ne-1]*uc[i, j, ne-1] - m[ne]*uc[i, j, ne]
        for i=1:nh, j=1:nc
    )
    # Dynamics
    ExaModels.@add_con(
        core,
        c6,
        - Duc[i,j,s] + g[s-1]*uc[i,j,s-1] - (m[s]+g[s])*uc[i,j,s]
        for i=1:nh, j=1:nc, s=2:ne-1
    )

    return core
end

@inline COPSBenchmark.marine_args(::ExaModelsBackend, nh; T = Float64) = (nh,)

@inline COPSBenchmark.marine_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.marine_recipe(b; T = T, backend = backend),
        COPSBenchmark.marine_args(b, nh; T = T)...;
        kwargs...,
    )
