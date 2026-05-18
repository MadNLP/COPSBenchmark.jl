# Rocket Steering Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.steering_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    a = T(100)  # Magnitude of force.
    # Bounds on the control
    u_min, u_max = -T(pi)/T(2), T(pi)/T(2)
    xs = zeros(T, 4)
    xf = [T(NaN), T(5), T(45), T(0)]
    half = T(0.5)
    inv_nh = T(1) / T(nh)

    function gen_x0(k, i)
        if i == 1 || i == 4
            return T(0)
        elseif i == 2
            return T(5)*T(k)*inv_nh
        elseif i == 3
            return T(45)*T(k)*inv_nh
        else
            return T(0)
        end
    end

    core = ExaModels.ExaCore(T; backend= backend, concrete = Val(true))

    ExaModels.@add_var(core, u, 1:nh+1; lvar = u_min, uvar =  u_max, start=T(0))   # control
    ExaModels.@add_var(core, x, 1:nh+1, 1:4; start=[gen_x0(i, j) for i=1:nh+1, j=1:4])     # state
    ExaModels.@add_var(core, tf, 1; start=T(1))                 # final time

    ExaModels.@add_obj(core, tf[1])

    ExaModels.@add_con(core, c0, tf[1]; lcon = T(0), ucon= T(Inf))
    # Dynamics
    ExaModels.@add_con(core, c1, -x[i+1,1] + x[i,1] + half*(tf[1]*inv_nh)*(x[i,3] + x[i+1,3]) for i=1:nh)
    ExaModels.@add_con(core, c2, -x[i+1,2] + x[i,2] + half*(tf[1]*inv_nh)*(x[i,4] + x[i+1,4]) for i=1:nh)
    ExaModels.@add_con(core, c3, -x[i+1,3] + x[i,3] + half*(tf[1]*inv_nh)*(a*cos(u[i]) + a*cos(u[i+1])) for i=1:nh)
    ExaModels.@add_con(core, c4, -x[i+1,4] + x[i,4] + half*(tf[1]*inv_nh)*(a*sin(u[i]) + a*sin(u[i+1])) for i=1:nh)
    # Boundary conditions
    ExaModels.@add_con(core, c5, -x[1, j] + s for (j,s) in enumerate(xs))
    ExaModels.@add_con(core, c6, -x[nh+1, j] + f for (j,f) in zip(2:4, xf[2:4]))

    return ExaModels.ExaModel(core; kwargs...)
end
