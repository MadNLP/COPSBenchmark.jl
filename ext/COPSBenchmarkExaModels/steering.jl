# Rocket Steering Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.steering_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    a = 100.0  # Magnitude of force.
    # Bounds on the control
    u_min, u_max = -pi/2.0, pi/2.0
    xs = zeros(4)
    xf = [NaN, 5.0, 45.0, 0.0]

    function gen_x0(k, i)
        if i == 1 || i == 4
            return 0.0
        elseif i == 2
            return 5*k/nh
        elseif i == 3
            return 45.0*k/nh
        else
            return 0.0
        end
    end

    core = ExaModels.ExaCore(T; backend= backend)

    ExaModels.@var(core, u, 1:nh+1; lvar = u_min, uvar =  u_max, start=0.0)   # control
    ExaModels.@var(core, x, 1:nh+1, 1:4; start=[gen_x0(i, j) for i=1:nh+1, j=1:4])     # state
    ExaModels.@var(core, tf, 1; start=1.0, lvar=0.0)                 # final time

    ExaModels.@obj(core, tf[1])

    # Dynamics
    ExaModels.@con(core, c1, -x[i+1,1] + x[i,1] + 0.5*(tf[1] / nh)*(x[i,3] + x[i+1,3]) for i=1:nh)
    ExaModels.@con(core, c2, -x[i+1,2] + x[i,2] + 0.5*(tf[1] / nh)*(x[i,4] + x[i+1,4]) for i=1:nh)
    ExaModels.@con(core, c3, -x[i+1,3] + x[i,3] + 0.5*(tf[1] / nh)*(a*cos(u[i]) + a*cos(u[i+1])) for i=1:nh)
    ExaModels.@con(core, c4, -x[i+1,4] + x[i,4] + 0.5*(tf[1] / nh)*(a*sin(u[i]) + a*sin(u[i+1])) for i=1:nh)
    # Boundary conditions
    ExaModels.@con(core, c5, -x[1, j] + s for (j,s) in enumerate(xs))
    ExaModels.@con(core, c6, -x[nh+1, j] + f for (j,f) in zip(2:4, xf[2:4]))

    return ExaModels.ExaModel(core; kwargs...)
end


