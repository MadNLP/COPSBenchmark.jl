# Goddard Rocket Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.rocket_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    h_0 = T(1)
    v_0 = T(0)
    m_0 = T(1)
    g_0 = T(1)
    T_c = T(3.5)
    h_c = T(500)
    v_c = T(620)
    m_c = T(0.6)
    half = T(0.5)
    inv_nh = T(1) / T(nh)

    c = half*sqrt(g_0 * h_0)
    m_f = m_c * m_0
    D_c = half * v_c * (m_0 / g_0)
    T_max = T_c * m_0 * g_0

    core = ExaModels.ExaCore(T; backend= backend, minimize=false, concrete = Val(true))

    # Precompute generator scalars so the broadcast closure stays isbits
    # (Metal rejects Type{T} captures). And bind the thrust-variable name as
    # Th rather than T to avoid shadowing the Type parameter T inside the
    # function — that shadowing caused Julia to Box the captured T value.
    s_v_start  = T[T(i)*inv_nh*(T(1) - T(i)*inv_nh) for i=0:nh]
    s_m_start  = T[(m_f - m_0)*(T(i)*inv_nh) + m_0 for i=0:nh]
    s_Th_start = T_max*half
    ExaModels.@add_var(core, h, 0:nh; start=h_0, lvar = h_0)
    ExaModels.@add_var(core, v, 0:nh; start=s_v_start, lvar = v_0)
    ExaModels.@add_var(core, m, 0:nh; start=s_m_start, lvar = m_f, uvar = m_0)
    ExaModels.@add_var(core, Th, 0:nh; start=s_Th_start, lvar = v_0, uvar = T_max)
    ExaModels.@add_var(core, step, 1; start=inv_nh, lvar = v_0)

    ExaModels.@add_obj(core, h[nh])

    # Dynamics
    ExaModels.@add_con(core, c1, - h[i] + h[i-1] + half * step[1] * (v[i] + v[i-1]) for i=1:nh)
    ExaModels.@add_con(core, c2, - v[i] + v[i-1] + half * step[1] * ((Th[i] - D_c*v[i]^2*exp(-h_c*(h[i] - h_0))/h_0 - m[i] * g_0 * (h_0 / h[i])^2) / m[i] + (Th[i-1] - D_c*v[i-1]^2*exp(-h_c*(h[i-1] - h_0))/h_0 - m[i-1] * g_0 * (h_0 / h[i-1])^2) / m[i-1]) for i=1:nh)
    ExaModels.@add_con(core, c3, - m[i] + m[i-1] + half * step[1] * (-Th[i]/c + -Th[i-1]/c) for i=1:nh)

    # Boundary ExaModels.constraints
    ExaModels.@add_con(core, c4, h[0] - h_0)
    ExaModels.@add_con(core, c5, v[0] - v_0)
    ExaModels.@add_con(core, c6, m[0] - m_0)
    ExaModels.@add_con(core, c7, m[nh] - m_f)

    return ExaModels.ExaModel(core; kwargs...)
end

function COPSBenchmark.rocket_core()
    args = ExaModels.ArgTracer()
    core = ExaCore(minimize = false, concrete = Val(true))
    h_0 = 1.0
    v_0 = 0.0
    m_0 = 1.0
    g_0 = 1.0
    T_c = 3.5
    h_c = 500.0
    v_c = 620.0
    m_c = 0.6
    half = 0.5
    inv_nh = 1 / args.nh

    c = half * sqrt(g_0 * h_0)
    m_f = m_c * m_0
    D_c = half * v_c * (m_0 / g_0)
    T_max = T_c * m_0 * g_0

    s_v_start = [i * inv_nh * (1.0 - i * inv_nh) for i = 0:args.nh]
    s_m_start = [(m_f - m_0) * (i * inv_nh) + m_0 for i = 0:args.nh]
    s_Th_start = T_max * half
    core, h = add_var(core, 0:args.nh; start = h_0, lvar = h_0)
    core, v = add_var(core, 0:args.nh; start = s_v_start, lvar = v_0)
    core, m = add_var(core, 0:args.nh; start = s_m_start, lvar = m_f, uvar = m_0)
    core, Th = add_var(core, 0:args.nh; start = s_Th_start, lvar = v_0, uvar = T_max)
    core, step = add_var(core, 1; start = inv_nh, lvar = v_0)

    core, _ = add_obj(core, h[k] for k in args.nh:args.nh)

    # Dynamics
    core, _ = add_con(core, -h[i] + h[i-1] + half * step[1] * (v[i] + v[i-1]) for i = 1:args.nh)
    core, _ = add_con(core, -v[i] + v[i-1] + half * step[1] * ((Th[i] - D_c * v[i]^2 * exp(-h_c * (h[i] - h_0)) / h_0 - m[i] * g_0 * (h_0 / h[i])^2) / m[i] + (Th[i-1] - D_c * v[i-1]^2 * exp(-h_c * (h[i-1] - h_0)) / h_0 - m[i-1] * g_0 * (h_0 / h[i-1])^2) / m[i-1]) for i = 1:args.nh)
    core, _ = add_con(core, -m[i] + m[i-1] + half * step[1] * (-Th[i] / c + -Th[i-1] / c) for i = 1:args.nh)

    # Boundary constraints
    core, _ = add_con(core, h[0] - h_0 for _ in 1:1)
    core, _ = add_con(core, v[0] - v_0 for _ in 1:1)
    core, _ = add_con(core, m[0] - m_0 for _ in 1:1)
    core, _ = add_con(core, m[k] - m_f for k in args.nh:args.nh)
    core
end
