# Electrons on a Sphere Problem.
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.elec_model(::JuMPBackend, np)
    Random.seed!(2713)

    theta = (2pi) .* rand(np)
    phi = pi .* rand(np)

    model = Model()
    @variable(model, x[i=1:np], start=cos(theta[i])*sin(phi[i]))
    @variable(model, y[i=1:np], start=sin(theta[i])*sin(phi[i]))
    @variable(model, z[i=1:np], start=cos(phi[i]))

    @expression(
        model,
        potential[i=1:np, j=i+1:np],
        COPSBenchmark.elec_coulomb_potential(x[i], y[i], z[i], x[j], y[j], z[j])
    )
    @objective(model, Min, sum(potential[i, j] for i in 1:np-1, j in i+1:np))

    @constraint(model, [i=1:np], x[i]^2 + y[i]^2 + z[i]^2 == 1)
    return model
end
