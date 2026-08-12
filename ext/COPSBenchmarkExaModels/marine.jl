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

    core, nh, d = ExaModels.ExaCore(T; backend = backend, nargs = Val(2))

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

@inline function COPSBenchmark.marine_args(::ExaModelsBackend, nh; T = Float64)
    ne = 8
    nm = 21
    tau = collect(T, range(T(0), T(10), 21))
    tf = tau[nm]
    h = tf / T(nh)
    t = T[T(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(T[
        20000.0, 17000.0, 10000.0, 15000.0, 12000.0,  9000.0,  7000.0,  3000.0,
        12445.0, 15411.0, 13040.0, 13338.0, 13484.0,  8426.0,  6615.0,  4022.0,
        7705.0, 13074.0, 14623.0, 11976.0, 12453.0,  9272.0,  6891.0,  5020.0,
        4664.0,  8579.0, 12434.0, 12603.0, 11738.0,  9710.0,  6821.0,  5722.0,
        2977.0,  7053.0, 11219.0, 11340.0, 13665.0,  8534.0,  6242.0,  5695.0,
        1769.0,  5054.0, 10065.0, 11232.0, 12112.0,  9600.0,  6647.0,  7034.0,
        943.0,  3907.0,  9473.0, 10334.0, 11115.0,  8826.0,  6842.0,  7348.0,
        581.0,  2624.0,  7421.0, 10297.0, 12427.0,  8747.0,  7199.0,  7684.0,
        355.0,  1744.0,  5369.0,  7748.0, 10057.0,  8698.0,  6542.0,  7410.0,
        223.0,  1272.0,  4713.0,  6869.0,  9564.0,  8766.0,  6810.0,  6961.0,
        137.0,   821.0,  3451.0,  6050.0,  8671.0,  8291.0,  6827.0,  7525.0,
        87.0,   577.0,  2649.0,  5454.0,  8430.0,  7411.0,  6423.0,  8388.0,
        49.0,   337.0,  2058.0,  4115.0,  7435.0,  7627.0,  6268.0,  7189.0,
        32.0,   228.0,  1440.0,  3790.0,  6474.0,  6658.0,  5859.0,  7467.0,
        17.0,   168.0,  1178.0,  3087.0,  6524.0,  5880.0,  5562.0,  7144.0,
        11.0,    99.0,   919.0,  2596.0,  5360.0,  5762.0,  4480.0,  7256.0,
        7.0,    65.0,   647.0,  1873.0,  4556.0,  5058.0,  4944.0,  7538.0,
        4.0,    44.0,   509.0,  1571.0,  4009.0,  4527.0,  4233.0,  6649.0,
        2.0,    27.0,   345.0,  1227.0,  3677.0,  4229.0,  3805.0,  6378.0,
        1.0,    20.0,   231.0,   934.0,  3197.0,  3695.0,  3159.0,  6454.0,
        1.0,    12.0,   198.0,   707.0,  2562.0,  3163.0,  3232.0,  5566.0,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = z[1, s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start = [v0[i, s] for i=1:nh, s=1:ne]
    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    return (nh, (; v_start, itr))
end

@inline COPSBenchmark.marine_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.marine_recipe(b; T = T, backend = backend),
        COPSBenchmark.marine_args(b, nh; T = T)...;
        kwargs...,
    )
