# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.glider_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    x_0    = COPSBenchmark.glider_x_0
    y_0    = COPSBenchmark.glider_y_0
    y_f    = COPSBenchmark.glider_y_f
    vx_0   = COPSBenchmark.glider_vx_0
    vx_f   = COPSBenchmark.glider_vx_f
    vy_0   = COPSBenchmark.glider_vy_0
    vy_f   = COPSBenchmark.glider_vy_f
    u_c    = COPSBenchmark.glider_u_c
    r_0    = COPSBenchmark.glider_r_0
    m      = COPSBenchmark.glider_m
    g      = COPSBenchmark.glider_g
    cd0    = COPSBenchmark.glider_cd0
    cd1    = COPSBenchmark.glider_cd1
    S      = COPSBenchmark.glider_S
    rho    = COPSBenchmark.glider_rho
    cL_min = COPSBenchmark.glider_cL_min
    cL_max = COPSBenchmark.glider_cL_max
    cL0    = cL_max / 2

    c = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(c, t_f, 1; lvar = 0.0, start = 1.0)
    ExaModels.@add_variable(c, x, nh+1; lvar = zeros(nh+1), start = [x_0 + vx_0*(k/nh) for k in 0:nh])
    ExaModels.@add_variable(c, y, nh+1; start = [y_0 + (k/nh)*(y_f - y_0) for k in 0:nh])
    ExaModels.@add_variable(c, vx, nh+1; lvar = zeros(nh+1), start = fill(vx_0, nh+1))
    ExaModels.@add_variable(c, vy, nh+1; start = fill(vy_0, nh+1))
    ExaModels.@add_variable(c, cL, nh+1; lvar = fill(cL_min,nh+1), uvar = fill(cL_max, nh+1), start = fill(cL0, nh+1))

    ExaModels.@add_expression(c, r, (x[i]/r_0 - 2.5)^2 for i in 1:nh+1)
    ExaModels.@add_expression(c, u, u_c*(1 - r[i])*exp(-r[i]) for i in 1:nh+1)
    ExaModels.@add_expression(c, w, vy[i] - u[i] for i in 1:nh+1)
    ExaModels.@add_expression(c, v, sqrt(vx[i]^2 + w[i]^2) for i in 1:nh+1)
    ExaModels.@add_expression(c, D, 0.5*(cd0+cd1*cL[i]^2)*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expression(c, L, 0.5*cL[i]*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expression(c, vx_dot, (-L[i]*(w[i]/v[i]) - D[i]*(vx[i]/v[i]))/m for i in 1:nh+1)
    ExaModels.@add_expression(c, vy_dot, (L[i]*(vx[i]/v[i]) - D[i]*(w[i]/v[i]))/m - g for i in 1:nh+1)

    ExaModels.@add_objective(c, -x[nh+1])

    ExaModels.@add_constraint(c, c1, x[j] - (x[j-1] + 0.5 * t_f[1]/nh * (vx[j] + vx[j-1])) for j in 2:nh+1)
    ExaModels.@add_constraint(c, c2, y[j] - (y[j-1] + 0.5 * t_f[1]/nh * (vy[j] + vy[j-1])) for j in 2:nh+1)
    ExaModels.@add_constraint(c, c3, vx[j] - (vx[j-1] + 0.5 * t_f[1]/nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:nh+1)
    ExaModels.@add_constraint(c, c4, vy[j] - (vy[j-1] + 0.5 * t_f[1]/nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:nh+1)

    ExaModels.@add_constraint(c, c5, x[1] - x_0)
    ExaModels.@add_constraint(c, c6, y[1] - y_0)
    ExaModels.@add_constraint(c, c7, y[nh+1] - y_f)
    ExaModels.@add_constraint(c, c8, vx[1] - vx_0)
    ExaModels.@add_constraint(c, c9, vx[nh+1] - vx_f)
    ExaModels.@add_constraint(c, c10, vy[1] - vy_0)
    ExaModels.@add_constraint(c, c11, vy[nh+1] - vy_f)

    ExaModels.ExaModel(c; kwargs...)
end
