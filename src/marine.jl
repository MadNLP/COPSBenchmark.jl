# Marine Population Dynamics Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function marine_model(nh)
    nc = 1     # number of collocation points
    ne = 8     # number of differential equations
    nm = 21    # number of measurements

    rho = [0.5]                            # roots of k-th degree Legendre polynomial
    tau = collect(range(0.0, 10.0, 21))    # times at which observations made
    tf = tau[nm]                           # ODEs defined in [0,tf]
    h = tf / nh                            # uniform interval length
    t = [(i-1)*h for i in 1:nh+1]          # partition

    # itau[i] is the largest integer k with t[k] <= tau[i]
    itau = Int[min(nh, floor(tau[i]/h)+1) for i in 1:nm]

    # Observation
    z = [
        20000.0 17000.0 10000.0 15000.0 12000.0  9000.0  7000.0  3000.0;
        12445.0 15411.0 13040.0 13338.0 13484.0  8426.0  6615.0  4022.0;
        7705.0 13074.0 14623.0 11976.0 12453.0  9272.0  6891.0  5020.0;
        4664.0  8579.0 12434.0 12603.0 11738.0  9710.0  6821.0  5722.0;
        2977.0  7053.0 11219.0 11340.0 13665.0  8534.0  6242.0  5695.0;
        1769.0  5054.0 10065.0 11232.0 12112.0  9600.0  6647.0  7034.0;
        943.0  3907.0  9473.0 10334.0 11115.0  8826.0  6842.0  7348.0;
        581.0  2624.0  7421.0 10297.0 12427.0  8747.0  7199.0  7684.0;
        355.0  1744.0  5369.0  7748.0 10057.0  8698.0  6542.0  7410.0;
        223.0  1272.0  4713.0  6869.0  9564.0  8766.0  6810.0  6961.0;
        137.0   821.0  3451.0  6050.0  8671.0  8291.0  6827.0  7525.0;
        87.0   577.0  2649.0  5454.0  8430.0  7411.0  6423.0  8388.0;
        49.0   337.0  2058.0  4115.0  7435.0  7627.0  6268.0  7189.0;
        32.0   228.0  1440.0  3790.0  6474.0  6658.0  5859.0  7467.0;
        17.0   168.0  1178.0  3087.0  6524.0  5880.0  5562.0  7144.0;
        11.0    99.0   919.0  2596.0  5360.0  5762.0  4480.0  7256.0;
        7.0    65.0   647.0  1873.0  4556.0  5058.0  4944.0  7538.0;
        4.0    44.0   509.0  1571.0  4009.0  4527.0  4233.0  6649.0;
        2.0    27.0   345.0  1227.0  3677.0  4229.0  3805.0  6378.0;
        1.0    20.0   231.0   934.0  3197.0  3695.0  3159.0  6454.0;
        1.0    12.0   198.0   707.0  2562.0  3163.0  3232.0  5566.0;
    ]

    v0 = zeros(nh, ne)
    # Starting-value
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = z[1, s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    model = Model()
    # Growth rates
    @variable(model, g[i=1:ne-1] >= 0.0)
    # Mortality rates
    @variable(model, m[i=1:ne] >= 0.0)
    # The collocation approximation u is defined by the parameters v and w.
    # uc and Duc are, respectively, u and u' evaluated at the collocation points.
    @variable(model, v[i=1:nh, s=1:ne], start=v0[i, s])
    @variable(model, w[i=1:nh, j=1:nc, s=1:ne], start=0.0)
    @variable(model, uc[i=1:nh, j=1:nc, s=1:ne], start=v0[i, s])
    @variable(model, Duc[i=1:nh, j=1:nc, s=1:ne], start=0.0)

    # error
    @expression(
        model,
        error[j=1:nm, s=1:ne],
        v[itau[j],s] + sum(w[itau[j],k,s]*(tau[j]-t[itau[j]])^k/(factorial(k)*h^(k-1)) for k in 1:nc) - z[j,s]
    )

    # L2 error
    @objective(model, Min, sum(error[j, s]^2 for j in 1:nm, s in 1:ne))

    # Collocation model
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
    # Continuity
    @constraint(
        model,
        [i=1:nh-1, s=1:ne],
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) == v[i+1, s]
    )
    # Boundary conditions
    @constraint(
        model,
        [i=1:nh, j=1:nc],
        Duc[i, j, 1] == -(m[1]+g[1])*uc[i, j, 1],
    )
    @constraint(
        model,
        [i=1:nh, j=1:nc],
        Duc[i, j, ne] == g[ne-1]*uc[i, j, ne-1] - m[ne]*uc[i, j, ne],
    )
    # Dynamics
    @constraint(
        model,
        [i=1:nh, j=1:nc, s=2:ne-1],
        Duc[i,j,s] == g[s-1]*uc[i,j,s-1] - (m[s]+g[s])*uc[i,j,s],
    )

    return model
end

