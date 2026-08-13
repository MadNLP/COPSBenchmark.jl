# Methanol-to-Hydrocarbons Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# `v0` is computed elaborately and then overwritten wholesale with 0.001, so
# every start in this model is that scalar -- nothing size-dependent survives.
# Only the measurement table, which depends on the interval length, is data.
@inline function COPSBenchmark.methanol_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    ne = 3
    np = 5
    nc = 3

    rho = T[0.11270166537926, 0.5, 0.88729833462074]
    tf = T(1.122)                         # tau[nm]; ODEs defined in [0,tf]
    zero_T = T(0)
    two_T = T(2)
    bc = T[1, 0, 0]

    c, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.methanol_data, nh)

    h = tf / nh                           # uniform interval length

    ExaModels.@add_var(c, theta, np; lvar = zero_T, start = fill(T(1), np))
    ExaModels.@add_var(c, v, nh, ne; start = T(0.001))
    ExaModels.@add_var(c, w, nh, nc, ne; start = zero_T)
    ExaModels.@add_var(c, uc, nh, nc, ne; start = T(0.001))
    ExaModels.@add_var(c, Duc, nh, nc, ne; start = zero_T)

    ExaModels.@add_obj(c, (v[itau,s] + sum(w[itau,k,s]*(tau-t)^k/(T(factorial(k))*h^(k-1)) for k in 1:nc) - z)^2 for (j,s,itau,tau,z,t) in d.con1_matrix)

    ExaModels.@add_con(
        c,
        c1,
        uc[i, j, s] - v[i,s] - h*sum(w[i,k,s]*(rho^k/T(factorial(k))) for k in 1:nc) for i=1:nh, (j,rho) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )

    ExaModels.@add_con(
        c,
        c2,
        Duc[i, j, s] - sum(w[i,k,s]*(rho^(k-1)/T(factorial(k-1))) for k in 1:nc) for i=1:nh, (j,rho) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )

    ExaModels.@add_con(
        c,
        c3,
        v[1, s] - bc for (s, bc) in [(s, bc[s]) for s in 1:ne]

    )

    ExaModels.@add_con(
        c,
        c4,
        v[i, s] + sum(w[i, j, s]*h/T(factorial(j)) for j in 1:nc) - v[i+1, s] for i=1:nh-1, s=1:ne
    )

    ExaModels.@add_con(
        c,
        c5,
        Duc[i,j,1] + ((two_T*theta[2] - (theta[1]*uc[i,j,2])/((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                         theta[3] + theta[4])*uc[i,j,1]) for i=1:nh, j=1:nc
    )

    ExaModels.@add_con(
        c,
        c6,
        Duc[i,j,2] - ((theta[1]*uc[i,j,1]*(theta[2]*uc[i,j,1]-uc[i,j,2]))/ ((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                     theta[3]*uc[i,j,1]) for i=1:nh, j=1:nc
    )

    ExaModels.@add_con(
        c,
        c7,
        Duc[i,j,3] - ((theta[1]*uc[i,j,1]*(uc[i,j,2]+theta[5]*uc[i,j,1]))/ ((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                        theta[4]*uc[i,j,1]) for i=1:nh, j=1:nc
    )

    return c
end

@inline COPSBenchmark.methanol_args(::ExaModelsBackend, nh; T = Float64) = (nh,)

@inline COPSBenchmark.methanol_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.methanol_recipe(b; T = T, backend = backend),
        COPSBenchmark.methanol_args(b, nh; T = T)...;
        kwargs...,
    )
