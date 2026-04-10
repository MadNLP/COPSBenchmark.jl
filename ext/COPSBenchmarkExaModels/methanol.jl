# Methanol-to-Hydrocarbons Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.methanol_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    ne = COPSBenchmark.methanol_ne
    np = COPSBenchmark.methanol_np
    nc = COPSBenchmark.methanol_nc
    nm = COPSBenchmark.methanol_nm
    rho = COPSBenchmark.methanol_rho
    bc  = COPSBenchmark.methanol_bc
    tau = COPSBenchmark.methanol_tau
    (; h, t, itau, z, v0) = COPSBenchmark.methanol_data(nh; T)

    con1_matrix = [(j, s, itau[j], tau[j], z[j,s], t[itau[j]]) for j in 1:nm, s in 1:ne]

    c = ExaModels.ExaCore(T; backend = backend)
    ExaModels.@add_variable(c, theta, np; lvar = 0, start = fill(1, np))
    ExaModels.@add_variable(c, v, nh, ne; start = v0)
    ExaModels.@add_variable(c, w, nh, nc, ne; start = 0)
    ExaModels.@add_variable(c, uc, nh, nc, ne; start = [v0[i,s] for i=1:nh, j=1:nc, s=1:ne])
    ExaModels.@add_variable(c, Duc, nh, nc, ne; start = 0)

    ExaModels.@add_objective(c, (v[itau_j,s] + sum(w[itau_j,k,s]*(tau_j-t_j)^k/(factorial(k)*h^(k-1)) for k in 1:nc) - z_js)^2 for (j,s,itau_j,tau_j,z_js,t_j) in con1_matrix)

    ExaModels.@add_constraint(
        c, c1,
        uc[i, j, s] - v[i,s] - h*sum(w[i,k,s]*(rhov^k/factorial(k)) for k in 1:nc) for i=1:nh, (j,rhov) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )
    ExaModels.@add_constraint(
        c, c2,
        Duc[i, j, s] - sum(w[i,k,s]*(rhov^(k-1)/factorial(k-1)) for k in 1:nc) for i=1:nh, (j,rhov) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )
    ExaModels.@add_constraint(c, c3, v[1, s] - bc_s for (s, bc_s) in [(s, bc[s]) for s in 1:ne])
    ExaModels.@add_constraint(
        c, c4,
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) - v[i+1, s] for i=1:nh-1, s=1:ne
    )
    ExaModels.@add_constraint(c, c5, Duc[i,j,1] - COPSBenchmark.methanol_ode1(theta, uc[i,j,1], uc[i,j,2]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(c, c6, Duc[i,j,2] - COPSBenchmark.methanol_ode2(theta, uc[i,j,1], uc[i,j,2]) for i=1:nh, j=1:nc)
    ExaModels.@add_constraint(c, c7, Duc[i,j,3] - COPSBenchmark.methanol_ode3(theta, uc[i,j,1], uc[i,j,2]) for i=1:nh, j=1:nc)

    return ExaModels.ExaModel(c; kwargs...)
end
