# Goddard Rocket Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.rocket_model(nh, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)
    h_0 = 1.0
    v_0 = 0.0
    m_0 = 1.0
    g_0 = 1.0
    T_c = 3.5
    h_c = 500.0
    v_c = 620.0
    m_c = 0.6

    c = 0.5*sqrt(g_0 * h_0)
    m_f = m_c * m_0
    D_c = 0.5 * v_c * (m_0 / g_0)
    T_max = T_c * m_0 * g_0

    core = ExaModels.ExaCore(T; backend= backend, minimize=false)

    h = ExaModels.variable(core,0:nh; start=1.0, lvar = 1.0)
    v = ExaModels.variable(core,0:nh; start=(i/nh*(1.0 - i/nh) for i=0:nh), lvar = 0.0)
    m = ExaModels.variable(core,0:nh; start=((m_f - m_0)*(i/nh) + m_0 for i=0:nh), lvar = m_f, uvar = m_0)
    T = ExaModels.variable(core,0:nh; start=T_max/2.0, lvar = 0.0, uvar = T_max)
    step = ExaModels.variable(core; start=1/nh, lvar = 0.0)

    ExaModels.objective(core, h[nh])

    # Dynamics
    ExaModels.constraint(core, - h[i] + h[i-1] + 0.5 * step * (v[i] + v[i-1]) for i=1:nh)
    ExaModels.constraint(core, - v[i] + v[i-1] + 0.5 * step * ((T[i] - D_c*v[i]^2*exp(-h_c*(h[i] - h_0))/h_0 - m[i] * g_0 * (h_0 / h[i])^2) / m[i] + (T[i-1] - D_c*v[i-1]^2*exp(-h_c*(h[i-1] - h_0))/h_0 - m[i-1] * g_0 * (h_0 / h[i-1])^2) / m[i-1]) for i=1:nh)
    ExaModels.constraint(core, - m[i] + m[i-1] + 0.5 * step * (-T[i]/c + -T[i-1]/c) for i=1:nh)

    # Boundary ExaModels.constraints
    ExaModels.constraint(core, h[0] - h_0)
    ExaModels.constraint(core, v[0] - v_0)
    ExaModels.constraint(core, m[0] - m_0)
    ExaModels.constraint(core, m[nh] - m_f)

    return ExaModels.ExaModel(core; kwargs...)
end

