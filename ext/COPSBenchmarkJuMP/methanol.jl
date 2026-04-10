# Methanol-to-Hydrocarbons Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.methanol_model(::JuMPBackend, nh)
    ne = COPSBenchmark.methanol_ne
    np = COPSBenchmark.methanol_np
    nc = COPSBenchmark.methanol_nc
    nm = COPSBenchmark.methanol_nm
    rho = COPSBenchmark.methanol_rho
    bc  = COPSBenchmark.methanol_bc
    tau = COPSBenchmark.methanol_tau
    (; h, t, itau, z, v0) = COPSBenchmark.methanol_data(nh)

    model = Model()

    @variable(model, theta[1:np] >= 0.0, start=1.0)
    @variable(model, v[i=1:nh, s=1:ne], start=v0[i, s])
    @variable(model, w[1:nh, 1:nc, 1:ne], start=0.0)
    @variable(model, uc[i=1:nh, j=1:nc, s=1:ne], start=v0[i, s])
    @variable(model, Duc[i=1:nh, j=1:nc, s=1:ne], start=0.0)

    @expression(
        model,
        error[j=1:nm, s=1:ne],
        v[itau[j],s] + sum(w[itau[j],k,s]*(tau[j]-t[itau[j]])^k/(factorial(k)*h^(k-1)) for k in 1:nc) - z[j,s]
    )
    @objective(model, Min, sum(error[j, s]^2 for j in 1:nm, s in 1:ne))

    @constraint(
        model,
        [i=1:nh, j=1:nc, s=1:ne],
        uc[i, j, s] == v[i,s] + h*sum(w[i,k,s]*(rho[j]^k/factorial(k)) for k in 1:nc),
    )
    @constraint(
        model,
        [i=1:nh, j=1:nc, s=1:ne],
        Duc[i, j, s] == sum(w[i,k,s]*(rho[j]^(k-1)/factorial(k-1)) for k in 1:nc),
    )
    @constraint(model, [s=1:ne], v[1, s] == bc[s])
    @constraint(
        model,
        [i=1:nh-1, s=1:ne],
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) == v[i+1, s]
    )
    @constraint(model, [i=1:nh, j=1:nc], Duc[i,j,1] == COPSBenchmark.methanol_ode1(theta, uc[i,j,1], uc[i,j,2]))
    @constraint(model, [i=1:nh, j=1:nc], Duc[i,j,2] == COPSBenchmark.methanol_ode2(theta, uc[i,j,1], uc[i,j,2]))
    @constraint(model, [i=1:nh, j=1:nc], Duc[i,j,3] == COPSBenchmark.methanol_ode3(theta, uc[i,j,1], uc[i,j,2]))

    return model
end
