# Goddard Rocket Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# The two varying starting points are data, not structure: they are
# comprehensions over the size, which has no symbolic form.  They are computed
# in `rocket_args` and supplied as arguments, so the core holds no function of
# this package's.  `inv_nh` does have a symbolic form and stays here.
@inline function COPSBenchmark.rocket_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    h_0 = T(1)
    v_0 = T(0)
    m_0 = T(1)
    g_0 = T(1)
    T_c = T(3.5)
    h_c = T(500)
    v_c = T(620)
    m_c = T(0.6)
    half = T(0.5)

    c = half*sqrt(g_0 * h_0)
    m_f = m_c * m_0
    D_c = half * v_c * (m_0 / g_0)
    T_max = T_c * m_0 * g_0

    core, nh = ExaModels.ExaCore(
        T; backend = backend, minimize = false, nargs = Val(1),
    )
    d = ExaModels.ArgCall(COPSBenchmark.rocket_data, (Val(T), nh))

    inv_nh = one(T) / nh
    s_Th_start = T_max*half

    ExaModels.@add_var(core, h, 0:nh; start=h_0, lvar = h_0)
    ExaModels.@add_var(core, v, 0:nh; start=d.v_start, lvar = v_0)
    ExaModels.@add_var(core, m, 0:nh; start=d.m_start, lvar = m_f, uvar = m_0)
    ExaModels.@add_var(core, Th, 0:nh; start=s_Th_start, lvar = v_0, uvar = T_max)
    ExaModels.@add_var(core, step, 1; start=inv_nh, lvar = v_0)

    # Indexes the last point, so the index arrives as data.
    ExaModels.@add_obj(core, h[k] for k in nh:nh)

    # Dynamics
    ExaModels.@add_con(core, c1, - h[i] + h[i-1] + half * step[1] * (v[i] + v[i-1]) for i=1:nh)
    ExaModels.@add_con(core, c2, - v[i] + v[i-1] + half * step[1] * ((Th[i] - D_c*v[i]^2*exp(-h_c*(h[i] - h_0))/h_0 - m[i] * g_0 * (h_0 / h[i])^2) / m[i] + (Th[i-1] - D_c*v[i-1]^2*exp(-h_c*(h[i-1] - h_0))/h_0 - m[i-1] * g_0 * (h_0 / h[i-1])^2) / m[i-1]) for i=1:nh)
    ExaModels.@add_con(core, c3, - m[i] + m[i-1] + half * step[1] * (-Th[i]/c + -Th[i-1]/c) for i=1:nh)

    # Boundary ExaModels.constraints
    ExaModels.@add_con(core, c4, h[0] - h_0)
    ExaModels.@add_con(core, c5, v[0] - v_0)
    ExaModels.@add_con(core, c6, m[0] - m_0)
    ExaModels.@add_con(core, c7, m[k] - m_f for k in nh:nh)

    return core
end

@inline COPSBenchmark.rocket_args(::ExaModelsBackend, nh) = (nh,)

@inline COPSBenchmark.rocket_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.rocket_recipe(b; T = T, backend = backend),
        COPSBenchmark.rocket_args(b, nh)...;
        kwargs...,
    )
