# Catalytic Cracking of Gas Oil Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# `h`, the interval length, has a symbolic form and stays in the recipe.  The
# tables that depend on it -- the measurement index map `itau`, the collocation
# times, and the starting values built from them -- do not, and travel as data.
@inline function COPSBenchmark.gasoil_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    nc = 4        # number of collocation points
    ne = 2        # number of differential equations
    np = 3        # number of ODE parameters
    nm = 21       # number of measurements

    rho = T[0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    tau = T[0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.150, 0.175, 0.20, 0.225, 0.250, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85, 0.95]
    tf = tau[nm]
    zero_T = T(0)
    z1 = T[1.0000, 0.0000]

    core, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode2(COPSBenchmark.gasoil_data, Val(T), nh)

    h = tf / nh

    # ODE parameters
    ExaModels.@add_var(core, theta, 1:np; lvar = zero_T, start=zero_T)
    ExaModels.@add_var(core, v, 1:nh, 1:ne; start=d.v_start)
    ExaModels.@add_var(core, w, 1:nh, 1:nc, 1:ne; start=zero_T)
    ExaModels.@add_var(core, uc, 1:nh, 1:nc, 1:ne; start=d.uc_start)
    ExaModels.@add_var(core, Duc, 1:nh, 1:nc, 1:ne; start=zero_T)

    itr2 = [(j,rho[j]) for j=1:nc]

    # L2 error
    ExaModels.@add_obj(core, (v[itauj,s] + sum(w[itauj,k,s]*(tauj-tj)^k/(T(factorial(k))*h^(k-1)) for k in 1:nc) - zjs)^2 for (j,s,itauj, tauj, tj, zjs) in d.itr)

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

    # Boundary
    itr3 = [(s,z1[s]) for s=1:ne]
    ExaModels.@add_con(core, c3, - v[1, s] + z1s for (s,z1s) in itr3)
    # Continuity
    ExaModels.@add_con(
        core,
        c4,
        v[i, s] + sum(w[i, j, s]*h/T(factorial(j)) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    ExaModels.@add_con(
        core,
        c5,
        - Duc[i, j, 1] - (theta[1]+theta[3])*uc[i, j, 1]^2
        for i=1:nh, j=1:nc
    )
    ExaModels.@add_con(
        core,
        c6,
        - Duc[i, j, 2] + theta[1]*uc[i,j,1]^2 - theta[2]*uc[i,j,2]
        for i=1:nh, j=1:nc
    )
    return core
end

@inline COPSBenchmark.gasoil_args(::ExaModelsBackend, nh) = (nh,)

@inline COPSBenchmark.gasoil_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.gasoil_recipe(b; T = T, backend = backend),
        COPSBenchmark.gasoil_args(b, nh)...;
        kwargs...,
    )
