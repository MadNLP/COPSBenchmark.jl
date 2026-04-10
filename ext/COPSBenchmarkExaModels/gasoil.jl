# Catalytic Cracking of Gas Oil Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.gasoil_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    nc = COPSBenchmark.gasoil_nc
    ne = COPSBenchmark.gasoil_ne
    np = COPSBenchmark.gasoil_np
    nm = COPSBenchmark.gasoil_nm
    rho = COPSBenchmark.gasoil_rho
    tau = COPSBenchmark.gasoil_tau
    (; h, t, itau, z, v0) = COPSBenchmark.gasoil_data(nh; T)

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(core, theta, 1:np; lvar = 0.0, start=0.0)
    ExaModels.@add_variable(core, v, 1:nh, 1:ne; start=[v0[i, s] for i=1:nh, s=1:ne])
    ExaModels.@add_variable(core, w, 1:nh, 1:nc, 1:ne; start=0.0)
    ExaModels.@add_variable(core, uc, 1:nh, 1:nc, 1:ne; start=[v0[i, s] for i=1:nh, j=1:nc, s=1:ne])
    ExaModels.@add_variable(core, Duc, 1:nh, 1:nc, 1:ne; start=0.0)

    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    itr2 = [(j,rho[j]) for j=1:nc]

    ExaModels.@add_objective(core, (v[itauj,s] + sum(w[itauj,k,s]*(tauj-tj)^k/(factorial(k)*h^(k-1)) for k in 1:nc) - zjs)^2 for (j,s,itauj, tauj, tj, zjs) in itr)

    ExaModels.@add_constraint(
        core, c1,
        - uc[i, j, s] + v[i,s] + h*sum(w[i,k,s]*(rhoj^k/factorial(k)) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    ExaModels.@add_constraint(
        core, c2,
        - Duc[i, j, s] + sum(w[i,k,s]*(rhoj^(k-1)/factorial(k-1)) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    # Boundary (use z[1,:] as initial condition)
    itr3 = [(s, z[1, s]) for s=1:ne]
    ExaModels.@add_constraint(core, c3, - v[1, s] + z1s for (s,z1s) in itr3)
    ExaModels.@add_constraint(
        core, c4,
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    ExaModels.@add_constraint(core, c5, - Duc[i, j, 1] + COPSBenchmark.gasoil_ode1(theta, uc[i, j, 1]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(core, c6, - Duc[i, j, 2] + COPSBenchmark.gasoil_ode2(theta, uc[i, j, 1], uc[i, j, 2]) for i=1:nh, j=1:nc)

    return ExaModels.ExaModel(core; kwargs...)
end
