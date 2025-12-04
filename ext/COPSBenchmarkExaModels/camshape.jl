# Cam Shape Problem
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.camshape_model(n, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)

    R_v = 1.0         # design parameter related to the valve shape
    R_max = 2.0       # maximum allowed radius of the cam
    R_min = 1.0       # minimum allowed radius of the cam
    alpha = 1.5       # curvature limit parameter

    d_theta = 2*pi/(5*(n+1))   # angle between discretization points

    core = ExaModels.ExaCore(T; backend= backend, minimize=false)

    # radius of the cam at discretization points
    r= ExaModels.variable(core, 1:n; lvar =R_min, uvar = R_max, start=(R_min+R_max)/2.0)

    ExaModels.objective(core, (pi*R_v)/n * r[i] for i in 1:n)

    # Convexity
    ExaModels.constraint(
        core,
        - r[i-1]*r[i] - r[i]*r[i+1] + 2*r[i-1]*r[i+1]*cos(d_theta) for i=2:n-1; lcon = -Inf, ucon = 0.0,
            )
    ExaModels.constraint(
        core,
        - R_min*r[1] - r[1]*r[2] + 2*R_min*r[2]*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.constraint(
        core,
        - R_min^2 - R_min*r[1] + 2*R_min*r[1]*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.constraint(
        core,
        - r[n-1]*r[n] - r[n]*R_max + 2*r[n-1]*R_max*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.constraint(
        core,
        - 2*R_max*r[n] + 2*r[n]^2*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    # Curvature
    ExaModels.constraint(
        core,
        (r[i+1] - r[i]) for i=1:n-1; lcon = -alpha*d_theta, ucon = alpha*d_theta,
    )
    ExaModels.constraint(
        core,
        (r[1] - R_min); lcon = -alpha*d_theta, ucon = alpha*d_theta
    )
    ExaModels.constraint(
        core,
        (R_max - r[n]); lcon = -alpha*d_theta, ucon = alpha*d_theta
    )

    return ExaModels.ExaModel(core; kwargs...)
end

