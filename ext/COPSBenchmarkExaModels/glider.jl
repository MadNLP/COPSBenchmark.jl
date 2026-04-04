# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.glider_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
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
    cd0 = 0.034
    cd1 = 0.069662
    S = 14.0
    rho = 1.13
    cL_min = 0.0
    cL_max = 1.4
    cL0 = cL_max / 2

    c = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@var(c, t_f, 1; lvar = 0.0, start = 1.0)
    ExaModels.@var(c, x, nh+1; lvar = zeros(nh+1), start = [x_0 + vx_0*(k/nh) for k in 0:nh])
    ExaModels.@var(c, y, nh+1; start = [y_0 + (k/nh)*(y_f - y_0) for k in 0:nh])
    ExaModels.@var(c, vx, nh+1; lvar = zeros(nh+1), start = fill(vx_0, nh+1))
    ExaModels.@var(c, vy, nh+1; start = fill(vy_0, nh+1))
    ExaModels.@var(c, cL, nh+1; lvar = fill(cL_min,nh+1), uvar = fill(cL_max, nh+1), start = fill(cL0, nh+1))

    # Expressions (matching JuMP @expressions)
    ExaModels.@expr(c, r, (x[i]/r_0 - 2.5)^2 for i in 1:nh+1)
    ExaModels.@expr(c, u, u_c*(1 - r[i])*exp(-r[i]) for i in 1:nh+1)
    ExaModels.@expr(c, w, vy[i] - u[i] for i in 1:nh+1)
    ExaModels.@expr(c, v, sqrt(vx[i]^2 + w[i]^2) for i in 1:nh+1)
    ExaModels.@expr(c, D, 0.5*(cd0+cd1*cL[i]^2)*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@expr(c, L, 0.5*cL[i]*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@expr(c, vx_dot, (-L[i]*(w[i]/v[i]) - D[i]*(vx[i]/v[i]))/m for i in 1:nh+1)
    ExaModels.@expr(c, vy_dot, (L[i]*(vx[i]/v[i]) - D[i]*(w[i]/v[i]))/m - g for i in 1:nh+1)

    ExaModels.@obj(c, -x[nh+1])

    # Dynamics
    ExaModels.@con(c, c1, x[j] - (x[j-1] + 0.5 * t_f[1]/nh * (vx[j] + vx[j-1])) for j in 2:nh+1)
    ExaModels.@con(c, c2, y[j] - (y[j-1] + 0.5 * t_f[1]/nh * (vy[j] + vy[j-1])) for j in 2:nh+1)
    ExaModels.@con(c, c3, vx[j] - (vx[j-1] + 0.5 * t_f[1]/nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:nh+1)
    ExaModels.@con(c, c4, vy[j] - (vy[j-1] + 0.5 * t_f[1]/nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:nh+1)

    # Boundary constraints
    ExaModels.@con(c, c5, x[1] - x_0)
    ExaModels.@con(c, c6, y[1] - y_0)
    ExaModels.@con(c, c7, y[nh+1] - y_f)
    ExaModels.@con(c, c8, vx[1] - vx_0)
    ExaModels.@con(c, c9, vx[nh+1] - vx_f)
    ExaModels.@con(c, c10, vy[1] - vy_0)
    ExaModels.@con(c, c11, vy[nh+1] - vy_f)

    ExaModels.ExaModel(c; kwargs...)
end
