# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.glider_model(::JuMPBackend, nh)
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
    c0     = COPSBenchmark.glider_cd0
    c1     = COPSBenchmark.glider_cd1
    S      = COPSBenchmark.glider_S
    rho    = COPSBenchmark.glider_rho
    cL_min = COPSBenchmark.glider_cL_min
    cL_max = COPSBenchmark.glider_cL_max

    model = Model()

    @variables(model, begin
        0 <= t_f,                       (start=1.0)
        0.0 <= x[k=0:nh],               (start=x_0 + vx_0*(k/nh))
        y[k=0:nh],                      (start=y_0 + (k/nh)*(y_f - y_0))
        0.0 <= vx[k=0:nh],              (start=vx_0)
        vy[k=0:nh],                     (start=vy_0)
        cL_min <= cL[k=0:nh] <= cL_max, (start=cL_max/2.0)
    end)

    @objective(model, Min, -x[nh])

    @expressions(model, begin
        step,           t_f / nh
        r[i=0:nh],      (x[i]/r_0 - 2.5)^2
        u[i=0:nh],      u_c*(1 - r[i])*exp(-r[i])
        w[i=0:nh],      vy[i] - u[i]
        v[i=0:nh],      sqrt(vx[i]^2 + w[i]^2)
        D[i=0:nh],      0.5*(c0+c1*cL[i]^2)*rho*S*v[i]^2
        L[i=0:nh],      0.5*cL[i]*rho*S*v[i]^2
        vx_dot[i=0:nh], (-L[i]*(w[i]/v[i]) - D[i]*(vx[i]/v[i]))/m
        vy_dot[i=0:nh], (L[i]*(vx[i]/v[i]) - D[i]*(w[i]/v[i]))/m - g
    end)

    @constraints(model, begin
        x_eqn[j=1:nh],  x[j] == x[j-1] + 0.5 * step * (vx[j] + vx[j-1])
        y_eqn[j=1:nh],  y[j] == y[j-1] + 0.5 * step * (vy[j] + vy[j-1])
        vx_eqn[j=1:nh], vx[j] == vx[j-1] + 0.5 * step * (vx_dot[j] + vx_dot[j-1])
        vy_eqn[j=1:nh], vy[j] == vy[j-1] + 0.5 * step * (vy_dot[j] + vy_dot[j-1])
    end)
    @constraints(model, begin
        x[0] == x_0
        y[0] == y_0
        y[nh] == y_f
        vx[0] == vx_0
        vx[nh] == vx_f
        vy[0] == vy_0
        vy[nh] == vy_f
    end)

    return model
end
