# Cam Shape Problem
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# The structure, with the number of discretization points left open.  `d_theta`
# depends on that number, so it cannot be a constant here — it is a deferred
# expression, and resolves to an ordinary Float64 when the model is built.  The
# rows that index the last point are written as generators over `n:n` so the
# index arrives as data rather than as a number the structure has to know.
@inline function COPSBenchmark.camshape_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    R_v = 1.0         # design parameter related to the valve shape
    R_max = 2.0       # maximum allowed radius of the cam
    R_min = 1.0       # minimum allowed radius of the cam
    alpha = 1.5       # curvature limit parameter

    core, n = ExaModels.ExaCore(
        T; backend = backend, minimize = false, nargs = Val(1),
    )

    d_theta = 2*pi/(5*(n+1))   # angle between discretization points

    # radius of the cam at discretization points
    ExaModels.@add_var(core, r, 1:n; lvar = R_min, uvar = R_max, start = (R_min+R_max)/2.0)

    ExaModels.@add_obj(core, (pi*R_v)/n * r[i] for i in 1:n)

    # Convexity
    ExaModels.@add_con(
        core,
        c1,
        - r[i-1]*r[i] - r[i]*r[i+1] + 2*r[i-1]*r[i+1]*cos(d_theta) for i=2:n-1; lcon = -Inf, ucon = 0.0,
            )
    ExaModels.@add_con(
        core,
        c2,
        - R_min*r[1] - r[1]*r[2] + 2*R_min*r[2]*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_con(
        core,
        c3,
        - R_min^2 - R_min*r[1] + 2*R_min*r[1]*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_con(
        core,
        c4,
        - r[k-1]*r[k] - r[k]*R_max + 2*r[k-1]*R_max*cos(d_theta) for k in n:n; lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_con(
        core,
        c5,
        - 2*R_max*r[k] + 2*r[k]^2*cos(d_theta) for k in n:n; lcon = -Inf, ucon = 0.0
    )
    # Curvature
    ExaModels.@add_con(
        core,
        c6,
        (r[i+1] - r[i]) for i=1:n-1; lcon = -alpha*d_theta, ucon = alpha*d_theta,
    )
    ExaModels.@add_con(
        core,
        c7,
        (r[1] - R_min); lcon = -alpha*d_theta, ucon = alpha*d_theta
    )
    ExaModels.@add_con(
        core,
        c8,
        (R_max - r[k]) for k in n:n; lcon = -alpha*d_theta, ucon = alpha*d_theta
    )

    return core
end

@inline COPSBenchmark.camshape_args(::ExaModelsBackend, n) = (n,)

@inline COPSBenchmark.camshape_model(b::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.camshape_recipe(b; T = T, backend = backend),
        COPSBenchmark.camshape_args(b, n)...;
        kwargs...,
    )
