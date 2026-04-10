# Cam Shape Problem
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.camshape_model(::ExaModelsBackend, n; T = Float64, backend = nothing, kwargs...)
    R_v   = COPSBenchmark.camshape_R_v
    R_max = COPSBenchmark.camshape_R_max
    R_min = COPSBenchmark.camshape_R_min
    alpha = COPSBenchmark.camshape_alpha

    d_theta = COPSBenchmark.camshape_d_theta(n)

    core = ExaModels.ExaCore(T; backend = backend, minimize=false)

    ExaModels.@add_variable(core, r, 1:n; lvar = R_min, uvar = R_max, start=(R_min+R_max)/2.0)

    ExaModels.@add_objective(core, (pi*R_v)/n * r[i] for i in 1:n)

    ExaModels.@add_constraint(
        core, c1,
        COPSBenchmark.camshape_convexity(r[i-1], r[i], r[i+1], d_theta) for i=2:n-1; lcon = -Inf, ucon = 0.0,
    )
    ExaModels.@add_constraint(
        core, c2,
        COPSBenchmark.camshape_convexity(R_min, r[1], r[2], d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_constraint(
        core, c3,
        - R_min^2 - R_min*r[1] + 2*R_min*r[1]*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_constraint(
        core, c4,
        COPSBenchmark.camshape_convexity(r[n-1], r[n], R_max, d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_constraint(
        core, c5,
        - 2*R_max*r[n] + 2*r[n]^2*cos(d_theta); lcon = -Inf, ucon = 0.0
    )
    ExaModels.@add_constraint(
        core, c6,
        (r[i+1] - r[i]) for i=1:n-1; lcon = -alpha*d_theta, ucon = alpha*d_theta,
    )
    ExaModels.@add_constraint(
        core, c7,
        (r[1] - R_min); lcon = -alpha*d_theta, ucon = alpha*d_theta
    )
    ExaModels.@add_constraint(
        core, c8,
        (R_max - r[n]); lcon = -alpha*d_theta, ucon = alpha*d_theta
    )

    return ExaModels.ExaModel(core; kwargs...)
end
