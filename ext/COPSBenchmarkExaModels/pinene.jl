# Isomerization of Alpha-Pinene Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004
#
@inline function COPSBenchmark.pinene_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    nc = COPSBenchmark.pinene_nc
    ne = COPSBenchmark.pinene_ne
    np = COPSBenchmark.pinene_np
    nm = COPSBenchmark.pinene_nm
    rho = COPSBenchmark.pinene_rho
    bc  = COPSBenchmark.pinene_bc
    tau = COPSBenchmark.pinene_tau
    (; h, t, itau, z, v0) = COPSBenchmark.pinene_data(nh; T)

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(core, theta, 1:np; lvar = 0.0, start=0.0)
    ExaModels.@add_variable(core, v, 1:nh, 1:ne; start=[v0[i, s] for i=1:nh, s=1:ne])
    ExaModels.@add_variable(core, w, 1:nh, 1:nc, 1:ne; start=0.0)
    ExaModels.@add_variable(core, uc, 1:nh, 1:nc, 1:ne; start=[v0[i,s] for i=1:nh, j=1:nc, s=1:ne])
    ExaModels.@add_variable(core, Duc, 1:nh, 1:nc, 1:ne; start=0.0)

    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    ExaModels.@add_objective(core, (v[it,s] + sum(w[it,k,s]*(tj-ti)^k/(factorial(k)*h^(k-1)) for k in 1:nc) - zjs)^2 for (j, s, it, tj, ti, zjs) in itr)

    itr2 = [(j,rho[j]) for j=1:nc]
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
    ExaModels.@add_constraint(core, c3, -v[1, s] + bcs for (s,bcs) in enumerate(bc))
    ExaModels.@add_constraint(
        core, c4,
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    ExaModels.@add_constraint(core, c5, -Duc[i,j,1] + COPSBenchmark.pinene_ode1(theta, uc[i,j,1]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(core, c6, -Duc[i,j,2] + COPSBenchmark.pinene_ode2(theta, uc[i,j,1]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(core, c7, -Duc[i,j,3] + COPSBenchmark.pinene_ode3(theta, uc[i,j,1], uc[i,j,3], uc[i,j,5]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(core, c8, -Duc[i,j,4] + COPSBenchmark.pinene_ode4(theta, uc[i,j,3]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(core, c9, -Duc[i,j,5] + COPSBenchmark.pinene_ode5(theta, uc[i,j,3], uc[i,j,5]) for i=1:nh, j=1:nc)

    return ExaModels.ExaModel(core; kwargs...)
end
