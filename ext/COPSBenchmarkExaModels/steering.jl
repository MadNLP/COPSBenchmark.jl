# Rocket Steering Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# The state's starting point is a matrix built by a local helper over the
# discretization, so it is data and moves to `steering_args`.  `u_min`/`u_max`
# and the boundary values do not depend on the size and stay here.
@inline function COPSBenchmark.steering_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    a = T(100)  # Magnitude of force.
    # Bounds on the control
    u_min, u_max = -T(pi)/T(2), T(pi)/T(2)
    xs = zeros(T, 4)
    xf = [T(NaN), T(5), T(45), T(0)]
    half = T(0.5)

    core, nh, x_start = ExaModels.ExaCore(T; backend = backend, nargs = Val(2))

    inv_nh = one(T) / nh

    ExaModels.@add_var(core, u, 1:nh+1; lvar = u_min, uvar =  u_max, start=T(0))   # control
    ExaModels.@add_var(core, x, 1:nh+1, 1:4; start=x_start)     # state
    ExaModels.@add_var(core, tf, 1; start=T(1))                 # final time

    ExaModels.@add_obj(core, tf[1])

    ExaModels.@add_con(core, c0, tf[1]; lcon = T(0), ucon= T(Inf))
    # Dynamics
    ExaModels.@add_con(core, c1, -x[i+1,1] + x[i,1] + half*(tf[1]*inv_nh)*(x[i,3] + x[i+1,3]) for i=1:nh)
    ExaModels.@add_con(core, c2, -x[i+1,2] + x[i,2] + half*(tf[1]*inv_nh)*(x[i,4] + x[i+1,4]) for i=1:nh)
    ExaModels.@add_con(core, c3, -x[i+1,3] + x[i,3] + half*(tf[1]*inv_nh)*(a*cos(u[i]) + a*cos(u[i+1])) for i=1:nh)
    ExaModels.@add_con(core, c4, -x[i+1,4] + x[i,4] + half*(tf[1]*inv_nh)*(a*sin(u[i]) + a*sin(u[i+1])) for i=1:nh)
    # Boundary conditions.  The final-row index depends on the size, so it comes
    # from a one-element set rather than being written into the structure.
    ExaModels.@add_con(core, c5, -x[1, j] + s for (j,s) in enumerate(xs))
    ExaModels.@add_con(core, c6, -x[k, j] + f for k in (nh+1):(nh+1), (j,f) in zip(2:4, xf[2:4]))

    return core
end

@inline function COPSBenchmark.steering_args(::ExaModelsBackend, nh; T = Float64)
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
    return (nh, [gen_x0(i, j) for i=1:nh+1, j=1:4])
end

@inline COPSBenchmark.steering_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.steering_recipe(b; T = T, backend = backend),
        COPSBenchmark.steering_args(b, nh; T = T)...;
        kwargs...,
    )
