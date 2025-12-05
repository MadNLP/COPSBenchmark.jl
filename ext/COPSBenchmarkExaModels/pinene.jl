# Isomerization of Alpha-Pinene Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004
#
function COPSBenchmark.pinene_model(nh, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)
    nc = 3        # number of collocation points
    ne = 5        # number of differential equations
    np = 5        # number of ODE parameters
    nm = 8        # number of measurements

    # roots of k-th degree Legendre polynomial
    rho = [0.11270166537926, 0.5, 0.88729833462074]
    # boundary conditions
    bc = [100.0, 0.0, 0.0, 0.0, 0.0]
    # times at which observations made
    tau = [1230.0, 3060.0, 4920.0, 7800.0, 10680.0, 15030.0, 22620.0, 36420.0]
    tf = tau[nm]                       # ODEs defined in [0,tf]
    h = tf / nh                        # uniform interval length
    t = [(i-1)*h for i in 1:nh+1]      # partition

    # itau[i] is the largest integer k with t[k] <= tau[i]
    itau = Int[min(nh, floor(tau[i]/h)+1) for i in 1:nm]

    # Observations
    z = [
        88.35    7.3     2.3     0.4     1.75;
        76.4    15.6     4.5     0.7     2.8;
        65.1    23.1     5.3     1.1     5.8;
        50.4    32.9     6.0     1.5     9.3;
        37.5    42.7     6.0     1.9    12.0;
        25.9    49.1     5.9     2.2    17.0;
        14.0    57.4     5.1     2.6    21.0;
        4.5    63.1     3.8     2.9    25.7;
    ]

    v0 = zeros(nh, ne)
    # Starting-value
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = bc[s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    core = ExaModels.ExaCore(T; backend= backend)

    theta = ExaModels.variable(core, 1:np; lvar = 0.0, start=0.0)
    # The collocation approximation u is defined by the parameters v and w.
    # uc and Duc are, respectively, u and u' evaluated at the collocation points.
    v = ExaModels.variable(core, 1:nh, 1:ne; start=[v0[i, s] for i=1:nh, s=1:ne])
    w = ExaModels.variable(core, 1:nh, 1:nc, 1:ne; start=0.0)
    uc= ExaModels.variable(core, 1:nh, 1:nc, 1:ne; start=[v0[i,s] for i=1:nh, j=1:nc, s=1:ne])
    Duc = ExaModels.variable(core, 1:nh, 1:nc, 1:ne; start=0.0)

    itr = [(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    # l2 error
    ExaModels.objective(core, (v[it,s] + sum(w[it,k,s]*(tj-ti)^k/(factorial(k)*h^(k-1)) for k in 1:nc) - zjs)^2 for (j, s, it, tj, ti, zjs) in itr)

    # Collocation model
    itr2 = [(j,rho[j]) for j=1:nc]
    ExaModels.constraint(
        core,
        - uc[i, j, s] + v[i,s] + h*sum(w[i,k,s]*(rhoj^k/factorial(k)) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    ExaModels.constraint(
        core,
        - Duc[i, j, s] + sum(w[i,k,s]*(rhoj^(k-1)/factorial(k-1)) for k in 1:nc)
        for i=1:nh, (j,rhoj) in itr2, s=1:ne
    )
    # Boundary
    ExaModels.constraint(core, -v[1, s] + bcs for (s,bcs) in enumerate(bc))
    # Continuity
    ExaModels.constraint(
        core,
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) - v[i+1, s]
        for i=1:nh-1, s=1:ne
    )
    ExaModels.constraint(core, -Duc[i,j,1] - (theta[1]+theta[2])*uc[i,j,1] for i=1:nh, j=1:nc)
    ExaModels.constraint(core, -Duc[i,j,2] + theta[1]*uc[i,j,1] for i=1:nh, j=1:nc)
    ExaModels.constraint(core, -Duc[i,j,3] + theta[2]*uc[i,j,1] - (theta[3]+theta[4])*uc[i,j,3] + theta[5]*uc[i,j,5] for i=1:nh, j=1:nc)
    ExaModels.constraint(core, -Duc[i,j,4] + theta[3]*uc[i,j,3] for i=1:nh, j=1:nc)
    ExaModels.constraint(core, -Duc[i,j,5] + theta[4]*uc[i,j,3] - theta[5]*uc[i,j,5] for i=1:nh, j=1:nc)

    return ExaModels.ExaModel(core; kwargs...)
end


