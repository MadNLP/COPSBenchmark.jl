# Isomerization of Alpha-Pinene Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004
#
# Same collocation template as gasoil/marine.  `h` keeps its symbolic form; the
# measurement index map and the starting values built from it are data.
@inline function COPSBenchmark.pinene_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    nc = 3        # number of collocation points
    ne = 5        # number of differential equations
    np = 5        # number of ODE parameters
    nm = 8        # number of measurements

    rho = T[0.11270166537926, 0.5, 0.88729833462074]
    bc = T[100, 0, 0, 0, 0]
    tf = T(36420)                           # tau[nm]
    zero_T = T(0)

    core, nh, d = ExaModels.ExaCore(T; backend = backend, nargs = Val(2))

    h = tf / nh                             # uniform interval length

    ExaModels.@add_var(core, theta, 1:np; lvar = zero_T, start=zero_T)
    # The collocation approximation u is defined by the parameters v and w.
    # uc and Duc are, respectively, u and u' evaluated at the collocation points.
    ExaModels.@add_var(core, v, 1:nh, 1:ne; start=d.v_start)
    ExaModels.@add_var(core, w, 1:nh, 1:nc, 1:ne; start=zero_T)
    ExaModels.@add_var(core, uc, 1:nh, 1:nc, 1:ne; start=d.uc_start)
    ExaModels.@add_var(core, Duc, 1:nh, 1:nc, 1:ne; start=zero_T)

    # l2 error
    ExaModels.@add_obj(core, (v[it,s] + sum(w[it,k,s]*(tj-ti)^k/(T(factorial(k))*h^(k-1)) for k in 1:nc) - zjs)^2 for (j, s, it, tj, ti, zjs) in d.itr)

    # Collocation model
    itr2 = [(j,rho[j]) for j=1:nc]
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
    ExaModels.@add_con(core, c3, -v[1, s] + bcs for (s,bcs) in enumerate(bc))
    # Continuity
    ExaModels.@add_con(
        core,
        c4,
        v[i, s] + sum(w[i, j, s]*h/T(factorial(j)) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    ExaModels.@add_con(core, c5, -Duc[i,j,1] - (theta[1]+theta[2])*uc[i,j,1] for i=1:nh, j=1:nc)
    ExaModels.@add_con(core, c6, -Duc[i,j,2] + theta[1]*uc[i,j,1] for i=1:nh, j=1:nc)
    ExaModels.@add_con(core, c7, -Duc[i,j,3] + theta[2]*uc[i,j,1] - (theta[3]+theta[4])*uc[i,j,3] + theta[5]*uc[i,j,5] for i=1:nh, j=1:nc)
    ExaModels.@add_con(core, c8, -Duc[i,j,4] + theta[3]*uc[i,j,3] for i=1:nh, j=1:nc)
    ExaModels.@add_con(core, c9, -Duc[i,j,5] + theta[4]*uc[i,j,3] - theta[5]*uc[i,j,5] for i=1:nh, j=1:nc)

    return core
end

@inline function COPSBenchmark.pinene_args(::ExaModelsBackend, nh; T = Float64)
    nc = 3
    ne = 5
    nm = 8
    bc = T[100, 0, 0, 0, 0]
    tau = T[1230, 3060, 4920, 7800, 10680, 15030, 22620, 36420]
    tf = tau[nm]
    h = tf / T(nh)
    t = T[T(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(T[
        88.35,  7.3,  2.3,  0.4,  1.75,
        76.4,  15.6,  4.5,  0.7,  2.8,
        65.1,  23.1,  5.3,  1.1,  5.8,
        50.4,  32.9,  6.0,  1.5,  9.3,
        37.5,  42.7,  6.0,  1.9, 12.0,
        25.9,  49.1,  5.9,  2.2, 17.0,
        14.0,  57.4,  5.1,  2.6, 21.0,
        4.5,  63.1,  3.8,  2.9, 25.7,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = bc[s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start  = [v0[i, s] for i=1:nh, s=1:ne]
    uc_start = [v0[i, s] for i=1:nh, j=1:nc, s=1:ne]
    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    return (nh, (; v_start, uc_start, itr))
end

@inline COPSBenchmark.pinene_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.pinene_recipe(b; T = T, backend = backend),
        COPSBenchmark.pinene_args(b, nh; T = T)...;
        kwargs...,
    )
