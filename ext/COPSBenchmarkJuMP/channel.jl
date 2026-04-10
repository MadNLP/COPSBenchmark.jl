# Flow in a Channel Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.channel_model(::JuMPBackend, nh)
    nc  = COPSBenchmark.channel_nc
    nd  = COPSBenchmark.channel_nd
    R   = COPSBenchmark.channel_R
    bc  = COPSBenchmark.channel_bc
    rho = COPSBenchmark.channel_rho

    (; h, v0) = COPSBenchmark.channel_v0(nh)

    model = Model()

    @variable(model, v[i=1:nh, j=1:nd], start=v0[i, j])
    @variable(model, w[i=1:nh, j=1:nc], start=0.0)

    @variable(model, uc[i=1:nh, j=1:nc, s=1:nd], start=v0[i, s])
    @variable(model, Duc[i=1:nh, j=1:nc, s=1:nd], start=0.0)

    # Constant objective
    @objective(model, Min, 1.0)

    # Collocation model
    @constraint(
        model,
        [i=1:nh, j=1:nc, s=1:nd],
        uc[i, j, s] == v[i,s] + h*sum(w[i,k]*(rho[j]^k/factorial(k)) for k in 1:nc),
    )
    @constraint(
        model,
        [i=1:nh, j=1:nc, s=1:nd],
        Duc[i, j, s] == sum(v[i,k]*((rho[j]*h)^(k-s)/factorial(k-s)) for k in s:nd) +
                            h^(nd-s+1) * sum(w[i, k]*(rho[j]^(k+nd-s)/factorial(k+nd-s)) for k in 1:nc)
    )
    # Boundary
    @constraint(model, bc_1, v[1, 1] == bc[1, 1])
    @constraint(model, bc_2, v[1, 2] == bc[2, 1])
    @constraint(
        model,
        bc_3,
        sum(v[nh, k]*(h^(k-1)/factorial(k-1)) for k in 1:nd) +
        h^nd * sum(w[nh, k]/factorial(k+nd-1) for k in 1:nc) == bc[1, 2]
    )
    @constraint(
        model,
        bc_4,
        sum(v[nh, k]*(h^(k-2)/factorial(k-2)) for k in 2:nd) +
        h^(nd-1) * sum(w[nh, k]/factorial(k+nd-2) for k in 1:nc) == bc[2, 2]
    )
    @constraint(
        model,
        continuity[i=1:nh-1, s=1:nd],
        sum(v[i, k]*(h^(k-s)/factorial(k-s)) for k in s:nd)
        + h^(nd-s+1)* sum(w[i, k]/factorial(k+nd-s) for k in 1:nc) == v[i+1, s]
    )
    @constraint(
        model,
        collocation[i=1:nh, j=1:nc],
        sum(w[i, k] * (rho[j]^(k-1)/factorial(k-1)) for k in 1:nc) ==
        COPSBenchmark.channel_collocation_rhs(R, Duc[i, j, 1], Duc[i, j, 2], Duc[i, j, 3], Duc[i, j, 4])
    )

    return model
end

