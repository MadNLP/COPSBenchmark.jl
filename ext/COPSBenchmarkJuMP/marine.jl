# Marine Population Dynamics Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.marine_model(::JuMPBackend, nh)
    nc = COPSBenchmark.marine_nc
    ne = COPSBenchmark.marine_ne
    nm = COPSBenchmark.marine_nm
    rho = COPSBenchmark.marine_rho
    tau = COPSBenchmark.marine_tau
    (; h, t, itau, z, v0) = COPSBenchmark.marine_data(nh)

    model = Model()

    @variable(model, g[i=1:ne-1] >= 0.0)
    @variable(model, m[i=1:ne] >= 0.0)
    @variable(model, v[i=1:nh, s=1:ne], start=v0[i, s])
    @variable(model, w[i=1:nh, j=1:nc, s=1:ne], start=0.0)
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
    @constraint(
        model,
        [i=1:nh-1, s=1:ne],
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) == v[i+1, s]
    )
    @constraint(model, [i=1:nh, j=1:nc], Duc[i, j, 1] == -(m[1]+g[1])*uc[i, j, 1])
    @constraint(model, [i=1:nh, j=1:nc], Duc[i, j, ne] == g[ne-1]*uc[i, j, ne-1] - m[ne]*uc[i, j, ne])
    @constraint(model, [i=1:nh, j=1:nc, s=2:ne-1], Duc[i,j,s] == g[s-1]*uc[i,j,s-1] - (m[s]+g[s])*uc[i,j,s])

    return model
end
