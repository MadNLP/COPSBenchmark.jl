# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.glider_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    # Design parameters (T-typed)
    x_0 = T(0)
    y_0 = T(1000)
    y_f = T(900)
    vx_0 = T(13.23)
    vx_f = T(13.23)
    vy_0 = T(-1.288)
    vy_f = T(-1.288)
    u_c = T(2.5)
    r_0 = T(100)
    m = T(100)
    g = T(9.81)
    cd0 = T(0.034)
    cd1 = T(0.069662)
    S = T(14)
    rho = T(1.13)
    cL_min = T(0)
    cL_max = T(1.4)
    cL0 = cL_max / T(2)
    half = T(0.5)
    inv_nh = T(1) / T(nh)
    r_offset = T(2.5)
    one_T = T(1)

    c = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))

    ExaModels.@add_var(c, t_f, 1; lvar = T(0), start = T(1))
    ExaModels.@add_var(c, x, nh+1; lvar = zeros(T, nh+1), start = [x_0 + vx_0*(T(k)*inv_nh) for k in 0:nh])
    ExaModels.@add_var(c, y, nh+1; start = [y_0 + (T(k)*inv_nh)*(y_f - y_0) for k in 0:nh])
    ExaModels.@add_var(c, vx, nh+1; lvar = zeros(T, nh+1), start = fill(vx_0, nh+1))
    ExaModels.@add_var(c, vy, nh+1; start = fill(vy_0, nh+1))
    ExaModels.@add_var(c, cL, nh+1; lvar = fill(cL_min,nh+1), uvar = fill(cL_max, nh+1), start = fill(cL0, nh+1))

    # Expressions (matching JuMP @expressions)
    ExaModels.@add_expr(c, r, (x[i]/r_0 - r_offset)^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, u, u_c*(one_T - r[i])*exp(-r[i]) for i in 1:nh+1)
    ExaModels.@add_expr(c, w, vy[i] - u[i] for i in 1:nh+1)
    ExaModels.@add_expr(c, v, sqrt(vx[i]^2 + w[i]^2) for i in 1:nh+1)
    ExaModels.@add_expr(c, D, half*(cd0+cd1*cL[i]^2)*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, L, half*cL[i]*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, vx_dot, (-L[i]*(w[i]/v[i]) - D[i]*(vx[i]/v[i]))/m for i in 1:nh+1)
    ExaModels.@add_expr(c, vy_dot, (L[i]*(vx[i]/v[i]) - D[i]*(w[i]/v[i]))/m - g for i in 1:nh+1)

    ExaModels.@add_obj(c, -x[nh+1])

    # Dynamics
    ExaModels.@add_con(c, c1, x[j] - (x[j-1] + half * t_f[1]*inv_nh * (vx[j] + vx[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c2, y[j] - (y[j-1] + half * t_f[1]*inv_nh * (vy[j] + vy[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c3, vx[j] - (vx[j-1] + half * t_f[1]*inv_nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c4, vy[j] - (vy[j-1] + half * t_f[1]*inv_nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:nh+1)

    # Boundary constraints
    ExaModels.@add_con(c, c5, x[1] - x_0)
    ExaModels.@add_con(c, c6, y[1] - y_0)
    ExaModels.@add_con(c, c7, y[nh+1] - y_f)
    ExaModels.@add_con(c, c8, vx[1] - vx_0)
    ExaModels.@add_con(c, c9, vx[nh+1] - vx_f)
    ExaModels.@add_con(c, c10, vy[1] - vy_0)
    ExaModels.@add_con(c, c11, vy[nh+1] - vy_f)

    ExaModels.ExaModel(c; kwargs...)
end

function COPSBenchmark.glider_core()
    args = ExaModels.ArgTracer()
    c = ExaCore(concrete = Val(true))
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
    half = 0.5
    inv_nh = 1 / args.nh
    r_offset = 2.5
    one_T = 1.0

    c, t_f = add_var(c, 1; lvar = 0.0, start = 1.0)
    c, x = add_var(c, args.nh+1; lvar = 0.0, start = [x_0 + vx_0 * (k * inv_nh) for k in 0:args.nh])
    c, y = add_var(c, args.nh+1; start = [y_0 + (k * inv_nh) * (y_f - y_0) for k in 0:args.nh])
    c, vx = add_var(c, args.nh+1; lvar = 0.0, start = fill(vx_0, args.nh+1))
    c, vy = add_var(c, args.nh+1; start = fill(vy_0, args.nh+1))
    c, cL = add_var(c, args.nh+1; lvar = fill(cL_min, args.nh+1), uvar = fill(cL_max, args.nh+1), start = fill(cL0, args.nh+1))

    # Expressions (matching JuMP @expressions)
    c, r = add_expr(c, (x[i] / r_0 - r_offset)^2 for i in 1:args.nh+1)
    c, u = add_expr(c, u_c * (one_T - r[i]) * exp(-r[i]) for i in 1:args.nh+1)
    c, w = add_expr(c, vy[i] - u[i] for i in 1:args.nh+1)
    c, v = add_expr(c, sqrt(vx[i]^2 + w[i]^2) for i in 1:args.nh+1)
    c, D = add_expr(c, half * (cd0 + cd1 * cL[i]^2) * rho * S * v[i]^2 for i in 1:args.nh+1)
    c, L = add_expr(c, half * cL[i] * rho * S * v[i]^2 for i in 1:args.nh+1)
    c, vx_dot = add_expr(c, (-L[i] * (w[i] / v[i]) - D[i] * (vx[i] / v[i])) / m for i in 1:args.nh+1)
    c, vy_dot = add_expr(c, (L[i] * (vx[i] / v[i]) - D[i] * (w[i] / v[i])) / m - g for i in 1:args.nh+1)

    c, _ = add_obj(c, -x[k] for k in args.nh+1:args.nh+1)

    # Dynamics
    c, _ = add_con(c, x[j] - (x[j-1] + half * t_f[1] * inv_nh * (vx[j] + vx[j-1])) for j in 2:args.nh+1)
    c, _ = add_con(c, y[j] - (y[j-1] + half * t_f[1] * inv_nh * (vy[j] + vy[j-1])) for j in 2:args.nh+1)
    c, _ = add_con(c, vx[j] - (vx[j-1] + half * t_f[1] * inv_nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:args.nh+1)
    c, _ = add_con(c, vy[j] - (vy[j-1] + half * t_f[1] * inv_nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:args.nh+1)

    # Boundary constraints
    c, _ = add_con(c, x[1] - x_0 for _ in 1:1)
    c, _ = add_con(c, y[1] - y_0 for _ in 1:1)
    c, _ = add_con(c, y[k] - y_f for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, vx[1] - vx_0 for _ in 1:1)
    c, _ = add_con(c, vx[k] - vx_f for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, vy[1] - vy_0 for _ in 1:1)
    c, _ = add_con(c, vy[k] - vy_f for k in args.nh+1:args.nh+1)
    c
end
