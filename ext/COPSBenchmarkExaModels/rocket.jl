# Goddard Rocket Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.rocket_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    h_0   = COPSBenchmark.rocket_h_0
    v_0   = COPSBenchmark.rocket_v_0
    m_0   = COPSBenchmark.rocket_m_0
    g_0   = COPSBenchmark.rocket_g_0
    h_c   = COPSBenchmark.rocket_h_c
    c     = COPSBenchmark.rocket_c
    m_f   = COPSBenchmark.rocket_m_f
    D_c   = COPSBenchmark.rocket_D_c
    T_max = COPSBenchmark.rocket_T_max

    core = ExaModels.ExaCore(T; backend = backend, minimize=false)

    ExaModels.@add_variable(core, h, 0:nh; start=1.0, lvar = 1.0)
    ExaModels.@add_variable(core, v, 0:nh; start=(i/nh*(1.0 - i/nh) for i=0:nh), lvar = 0.0)
    ExaModels.@add_variable(core, m, 0:nh; start=((m_f - m_0)*(i/nh) + m_0 for i=0:nh), lvar = m_f, uvar = m_0)
    ExaModels.@add_variable(core, T, 0:nh; start=T_max/2.0, lvar = 0.0, uvar = T_max)
    ExaModels.@add_variable(core, step, 1; start=1/nh, lvar = 0.0)

    ExaModels.@add_objective(core, h[nh])

    ExaModels.@add_constraint(core, c1, - h[i] + h[i-1] + 0.5 * step[1] * (v[i] + v[i-1]) for i=1:nh)
    ExaModels.@add_constraint(core, c2, - v[i] + v[i-1] + 0.5 * step[1] * ((T[i] - D_c*v[i]^2*exp(-h_c*(h[i] - h_0))/h_0 - m[i] * g_0 * (h_0 / h[i])^2) / m[i] + (T[i-1] - D_c*v[i-1]^2*exp(-h_c*(h[i-1] - h_0))/h_0 - m[i-1] * g_0 * (h_0 / h[i-1])^2) / m[i-1]) for i=1:nh)
    ExaModels.@add_constraint(core, c3, - m[i] + m[i-1] + 0.5 * step[1] * (-T[i]/c + -T[i-1]/c) for i=1:nh)

    ExaModels.@add_constraint(core, c4, h[0] - h_0)
    ExaModels.@add_constraint(core, c5, v[0] - v_0)
    ExaModels.@add_constraint(core, c6, m[0] - m_0)
    ExaModels.@add_constraint(core, c7, m[nh] - m_f)

    return ExaModels.ExaModel(core; kwargs...)
end
