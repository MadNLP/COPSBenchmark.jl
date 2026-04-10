# Catalytic Cracking of Gas Oil Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.gasoil_model(::JuMPBackend, nh)
    nc = COPSBenchmark.gasoil_nc
    ne = COPSBenchmark.gasoil_ne
    np = COPSBenchmark.gasoil_np
    nm = COPSBenchmark.gasoil_nm
    rho = COPSBenchmark.gasoil_rho
    tau = COPSBenchmark.gasoil_tau
    (; h, t, itau, z, v0) = COPSBenchmark.gasoil_data(nh)

    model = Model()

    @variable(model, theta[1:np] >= 0.0, start=0.0)
    @variable(model, v[i=1:nh, s=1:ne], start=v0[i, s])
    @variable(model, w[1:nh, 1:nc, 1:ne], start=0.0)
    @variable(model, uc[i=1:nh, j=1:nc, s=1:ne], start=v0[i,s])
    @variable(model, Duc[i=1:nh, j=1:nc, s=1:ne], start=0.0)

    @expression(
        model,
        error[j=1:nm, s=1:ne],
        v[itau[j],s] + sum(w[itau[j],k,s]*(tau[j]-t[itau[j]])^k/(factorial(k)*h^(k-1)) for k in 1:nc) - z[j,s],
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
    @constraint(model, [s=1:ne], v[1, s] == z[1, s])
    @constraint(
        model,
        [i=1:nh-1, s=1:ne],
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) == v[i+1, s],
    )
    @constraint(model, [i=1:nh, j=1:nc], Duc[i, j, 1] == COPSBenchmark.gasoil_ode1(theta, uc[i, j, 1]))
    @constraint(model, [i=1:nh, j=1:nc], Duc[i, j, 2] == COPSBenchmark.gasoil_ode2(theta, uc[i, j, 1], uc[i, j, 2]))
    return model
end
