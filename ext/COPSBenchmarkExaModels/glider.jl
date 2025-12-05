# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.glider_model(nh, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)
    # Design parameters
    x_0 = 0.0
    y_0 = 1000.0
    y_f = 900.0
    vx_0 = 13.23
    vx_f = 13.23
    vy_0 = -1.288
    vy_f = -1.288
    u_c = 2.5
    r_0 = 100.0
    m = 100.0
    g = 9.81
    c0 = 0.034
    c1 = 0.069662
    S = 14.0
    rho = 1.13
    cL_min = 0.0
    cL_max = 1.4
    cL0 = cL_max / 2

    # Initial position
    x0 = [x_0 + vx_0*(k/nh) for k in 0:nh]
    r0 = zeros(nh+1)
    u0 = zeros(nh+1)
    w0 = [vy_0 - u0[i] for i in 1:nh+1]
    v0 = [sqrt(vx_0^2 + w0[i]^2) for i in 1:nh+1]
    D0 = [0.5*(c0+c1*cL0^2)*rho*S*v0[i]^2 for i in 1:nh+1]
    L0 = [0.5*cL0*rho*S*v0[i]^2 for i in 1:nh+1]
    vxdot0 = [-L0[i]*w0[i]/v0[i]/m - D0[i]*vx_0/v0[i]/m for i in 1:nh+1]
    vydot0 = [L0[i]*vx_0/v0[i]/m - D0[i]*w0[i]/v0[i]/m for i in 1:nh+1]

    c = ExaModels.ExaCore(T; backend = backend)

    t_f = ExaModels.variable(c, 1; lvar = 0.1, start = 1.0)
    x = ExaModels.variable(c, nh+1; lvar = zeros(nh+1), start = x0)
    y = ExaModels.variable(c, nh+1; start = [y_0 + (k/nh)*(y_f - y_0) for k in 0:nh])
    vx = ExaModels.variable(c, nh+1; lvar = zeros(nh+1), start = fill(vx_0, nh+1))
    vy = ExaModels.variable(c, nh+1; start = fill(vy_0, nh+1))
    cL = ExaModels.variable(c, nh+1; lvar = fill(cL_min,nh+1), uvar = fill(cL_max, nh+1), start = fill(cL0, nh+1))

    # Functions that define the glider (lifted variable space)
    r = ExaModels.variable(c, nh+1, start=r0)
    u = ExaModels.variable(c, nh+1, start=u0)
    w = ExaModels.variable(c, nh+1, start=w0)
    v = ExaModels.variable(c, nh+1, start=v0)
    D = ExaModels.variable(c, nh+1, start=D0)
    L = ExaModels.variable(c, nh+1, start=L0)
    vx_dot = ExaModels.variable(c, nh+1, start=vxdot0)
    vy_dot = ExaModels.variable(c, nh+1, start=vydot0)

    ExaModels.objective(c, -x[nh+1])

    # Equations of motion.
    ExaModels.constraint(c, r[i] - (x[i]/r_0 - 2.5)^2 for i in 1:nh+1)
    ExaModels.constraint(c, u[i] - u_c*(1 - r[i])*exp(-r[i]) for i in 1:nh+1)
    ExaModels.constraint(c, w[i] - vy[i] + u[i] for i in 1:nh+1)
    ExaModels.constraint(c, v[i] - sqrt(vx[i]^2 + w[i]^2) for i in 1:nh+1)
    ExaModels.constraint(c, D[i] - 0.5*(c0+c1*cL[i]^2)*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.constraint(c, L[i] - 0.5*cL[i]*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.constraint(c, vx_dot[i] + L[i]*w[i]/v[i]/m + D[i]*vx[i]/v[i]/m for i in 1:nh+1)
    ExaModels.constraint(c, vy_dot[i] - L[i]*vx[i]/v[i]/m + D[i]*w[i]/v[i]/m + g for i in 1:nh+1)

    # Collocations
    ExaModels.constraint(
        c,
        x[j] - (x[j-1] + 0.5 * t_f[1]/nh * (vx[j] + vx[j-1])) for j in 2:nh+1
    )
    ExaModels.constraint(
        c,
        y[j] - (y[j-1] + 0.5 * t_f[1]/nh * (vy[j] + vy[j-1])) for j in 2:nh+1
    )
    ExaModels.constraint(
        c,
        vx[j] - (vx[j-1] + 0.5 * t_f[1]/nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:nh+1
    )
    ExaModels.constraint(
        c,
        vy[j] - (vy[j-1] + 0.5 * t_f[1]/nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:nh+1
    )

    ExaModels.constraint(
        c,
        x[1] - x_0
    )

    ExaModels.constraint(
        c,
        y[1] - y_0
    )

    ExaModels.constraint(
        c,
        y[nh+1] - y_f
    )

    ExaModels.constraint(
        c,
        vx[1] - vx_0
    )

    ExaModels.constraint(
        c,
        vx[nh+1] - vx_f
    )

    ExaModels.constraint(
        c,
        vy[1] - vy_0
    )

    ExaModels.constraint(
        c,
        vy[nh+1] - vy_f
    )

    ExaModels.ExaModel(c; kwargs...)
end

