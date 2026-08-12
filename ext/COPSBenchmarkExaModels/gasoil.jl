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

    core, nh, d = ExaModels.ExaCore(T; backend = backend, nargs = Val(2))

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

@inline function COPSBenchmark.gasoil_args(::ExaModelsBackend, nh; T = Float64)
    nc = 4
    ne = 2
    nm = 21
    rho = T[0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    bc = [1, 1, 2, 0]
    tau = T[0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.150, 0.175, 0.20, 0.225, 0.250, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85, 0.95]
    tf = tau[nm]
    h = tf / T(nh)
    t = T[T(i-1)*h for i in 1:nh+1]

    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(T[
        1.0000, 0.0000,
        0.8105, 0.2000,
        0.6208, 0.2886,
        0.5258, 0.3010,
        0.4345, 0.3215,
        0.3903, 0.3123,
        0.3342, 0.2716,
        0.3034, 0.2551,
        0.2735, 0.2258,
        0.2405, 0.1959,
        0.2283, 0.1789,
        0.2071, 0.1457,
        0.1669, 0.1198,
        0.1530, 0.0909,
        0.1339, 0.0719,
        0.1265, 0.0561,
        0.1200, 0.0460,
        0.0990, 0.0280,
        0.0870, 0.0190,
        0.0770, 0.0140,
        0.0690, 0.0100,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = T(bc[s])
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start  = [v0[i, s] for i = 1:nh, s = 1:ne]
    uc_start = [v0[i, s] for i = 1:nh, j = 1:nc, s = 1:ne]
    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]

    return (nh, (; v_start, uc_start, itr))
end

@inline COPSBenchmark.gasoil_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.gasoil_recipe(b; T = T, backend = backend),
        COPSBenchmark.gasoil_args(b, nh; T = T)...;
        kwargs...,
    )
