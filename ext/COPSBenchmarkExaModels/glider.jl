# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.glider_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
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
    r_offset = T(2.5)
    one_T = T(1)

    c, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgCall(COPSBenchmark.glider_data, (Val(T), nh))

    inv_nh = one(T) / nh

    ExaModels.@add_var(c, t_f, 1; lvar = T(0), start = T(1))
    ExaModels.@add_var(c, x, nh+1; lvar = T(0), start = d.x_start)
    ExaModels.@add_var(c, y, nh+1; start = d.y_start)
    ExaModels.@add_var(c, vx, nh+1; lvar = T(0), start = vx_0)
    ExaModels.@add_var(c, vy, nh+1; start = vy_0)
    ExaModels.@add_var(c, cL, nh+1; lvar = cL_min, uvar = cL_max, start = cL0)

    # Expressions (matching JuMP @expressions)
    ExaModels.@add_expr(c, r, (x[i]/r_0 - r_offset)^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, u, u_c*(one_T - r[i])*exp(-r[i]) for i in 1:nh+1)
    ExaModels.@add_expr(c, w, vy[i] - u[i] for i in 1:nh+1)
    ExaModels.@add_expr(c, v, sqrt(vx[i]^2 + w[i]^2) for i in 1:nh+1)
    ExaModels.@add_expr(c, D, half*(cd0+cd1*cL[i]^2)*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, L, half*cL[i]*rho*S*v[i]^2 for i in 1:nh+1)
    ExaModels.@add_expr(c, vx_dot, (-L[i]*(w[i]/v[i]) - D[i]*(vx[i]/v[i]))/m for i in 1:nh+1)
    ExaModels.@add_expr(c, vy_dot, (L[i]*(vx[i]/v[i]) - D[i]*(w[i]/v[i]))/m - g for i in 1:nh+1)

    ExaModels.@add_obj(c, -x[k] for k in (nh+1):(nh+1))

    # Dynamics
    ExaModels.@add_con(c, c1, x[j] - (x[j-1] + half * t_f[1]*inv_nh * (vx[j] + vx[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c2, y[j] - (y[j-1] + half * t_f[1]*inv_nh * (vy[j] + vy[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c3, vx[j] - (vx[j-1] + half * t_f[1]*inv_nh * (vx_dot[j] + vx_dot[j-1])) for j in 2:nh+1)
    ExaModels.@add_con(c, c4, vy[j] - (vy[j-1] + half * t_f[1]*inv_nh * (vy_dot[j] + vy_dot[j-1])) for j in 2:nh+1)

    # Boundary constraints
    ExaModels.@add_con(c, c5, x[1] - x_0)
    ExaModels.@add_con(c, c6, y[1] - y_0)
    ExaModels.@add_con(c, c7, y[k] - y_f for k in (nh+1):(nh+1))
    ExaModels.@add_con(c, c8, vx[1] - vx_0)
    ExaModels.@add_con(c, c9, vx[k] - vx_f for k in (nh+1):(nh+1))
    ExaModels.@add_con(c, c10, vy[1] - vy_0)
    ExaModels.@add_con(c, c11, vy[k] - vy_f for k in (nh+1):(nh+1))

    return c
end

@inline COPSBenchmark.glider_args(::ExaModelsBackend, nh) = (nh,)

@inline COPSBenchmark.glider_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.glider_recipe(b; T = T, backend = backend),
        COPSBenchmark.glider_args(b, nh)...;
        kwargs...,
    )
