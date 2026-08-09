# Robot Arm Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.robot_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)

    # total length of arm
    L = T(5)

    # Upper bounds on the controls
    max_u_rho = T(1)
    max_u_the = T(1)
    max_u_phi = T(1)
    # Initial positions of the length and the angles for the robot arm
    rho0 = T(4.5)
    pi_T = T(pi)
    phi0 = pi_T / T(4)
    half = T(0.5)
    third = T(1) / T(3)
    inv_nh = T(1) / T(nh)
    two_pi_3 = T(2) * pi_T / T(3)
    four_pi_3 = T(4) * pi_T / T(3)
    zero_T = T(0)

    core = ExaModels.ExaCore(T; backend= backend, concrete = Val(true))

    ExaModels.@add_var(core, rho, nh+1; start=rho0, lvar = zero_T, uvar = L)
    ExaModels.@add_var(core, the, nh+1; start=[two_pi_3*(T(k)*inv_nh)^2 for k=1:nh+1], lvar = -pi_T, uvar = pi_T)
    ExaModels.@add_var(core, phi, nh+1; start=phi0, lvar = zero_T, uvar = pi_T)
    # Derivatives
    ExaModels.@add_var(core, rho_dot, nh+1; start=zero_T)
    ExaModels.@add_var(core, the_dot, nh+1; start=[four_pi_3*(T(k)*inv_nh) for k=1:nh+1])
    ExaModels.@add_var(core, phi_dot, nh+1; start=zero_T)
    # Control
    ExaModels.@add_var(core, u_rho, nh+1; start=zero_T, lvar = -max_u_rho, uvar = max_u_rho)
    ExaModels.@add_var(core, u_the, nh+1; start=zero_T, lvar = -max_u_the, uvar = max_u_the)
    ExaModels.@add_var(core, u_phi, nh+1; start=zero_T, lvar = -max_u_phi, uvar = max_u_phi)
    # Final time
    ExaModels.@add_var(core, tf, 1; start = T(1), lvar = zero_T)

    # Expressions (matching JuMP @expressions)
    ExaModels.@add_expr(core, I_the, ((L-rho[i])^3+rho[i]^3)*(sin(phi[i]))^2*third for i=1:nh+1)
    ExaModels.@add_expr(core, I_phi, ((L-rho[i])^3+rho[i]^3)*third for i=1:nh+1)

    ExaModels.@add_obj(core, tf[1] for i=1:1)

    # Dynamics
    ExaModels.@add_con(core, c1, - rho[j] + rho[j-1] + half * tf[1]*inv_nh * (rho_dot[j] + rho_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_con(core, c2, - phi[j] + phi[j-1] + half * tf[1]*inv_nh * (phi_dot[j] + phi_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_con(core, c3, - the[j] + the[j-1] + half * tf[1]*inv_nh * (the_dot[j] + the_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_con(core, c4, - rho_dot[j] + rho_dot[j-1] + half * tf[1]*inv_nh * (u_rho[j] + u_rho[j-1]) / L for j=2:nh+1)
    ExaModels.@add_con(core, c5, - the_dot[j] + the_dot[j-1] + half * tf[1]*inv_nh * (u_the[j] / I_the[j] + u_the[j-1] / I_the[j-1]) for j=2:nh+1)
    ExaModels.@add_con(core, c6, - phi_dot[j] + phi_dot[j-1] + half * tf[1]*inv_nh * (u_phi[j] / I_phi[j] + u_phi[j-1] / I_phi[j-1]) for j=2:nh+1)

    # Boundary conditions
    ExaModels.@add_con(core, c7, - rho[1] + rho0)
    ExaModels.@add_con(core, c8, - the[1] + zero_T)
    ExaModels.@add_con(core, c9, - phi[1] + phi0)
    ExaModels.@add_con(core, c10, - rho[nh+1] + rho0)
    ExaModels.@add_con(core, c11, - the[nh+1] + two_pi_3)
    ExaModels.@add_con(core, c12, - phi[nh+1] + phi0)
    ExaModels.@add_con(core, c13, - rho_dot[1] + zero_T)
    ExaModels.@add_con(core, c14, - the_dot[1] + zero_T)
    ExaModels.@add_con(core, c15, - phi_dot[1] + zero_T)
    ExaModels.@add_con(core, c16, - rho_dot[nh+1] + zero_T)
    ExaModels.@add_con(core, c17, - the_dot[nh+1] + zero_T)
    ExaModels.@add_con(core, c18, - phi_dot[nh+1] + zero_T)

    return ExaModels.ExaModel(core; kwargs...)
end

function COPSBenchmark.robot_core()
    args = ExaModels.ArgTracer()
    c = ExaCore(concrete = Val(true))
    # total length of arm
    L = 5.0
    # Upper bounds on the controls
    max_u_rho = 1.0
    max_u_the = 1.0
    max_u_phi = 1.0
    # Initial positions of the length and the angles for the robot arm
    rho0 = 4.5
    pi_T = pi
    phi0 = pi_T / 4
    half = 0.5
    third = 1 / 3
    inv_nh = 1 / args.nh
    two_pi_3 = 2 * pi_T / 3
    four_pi_3 = 4 * pi_T / 3
    zero_T = 0.0

    c, rho = add_var(c, args.nh+1; start = rho0, lvar = zero_T, uvar = L)
    c, the = add_var(c, args.nh+1; start = [two_pi_3 * (k * inv_nh)^2 for k = 1:args.nh+1], lvar = -pi_T, uvar = pi_T)
    c, phi = add_var(c, args.nh+1; start = phi0, lvar = zero_T, uvar = pi_T)
    # Derivatives
    c, rho_dot = add_var(c, args.nh+1; start = zero_T)
    c, the_dot = add_var(c, args.nh+1; start = [four_pi_3 * (k * inv_nh) for k = 1:args.nh+1])
    c, phi_dot = add_var(c, args.nh+1; start = zero_T)
    # Control
    c, u_rho = add_var(c, args.nh+1; start = zero_T, lvar = -max_u_rho, uvar = max_u_rho)
    c, u_the = add_var(c, args.nh+1; start = zero_T, lvar = -max_u_the, uvar = max_u_the)
    c, u_phi = add_var(c, args.nh+1; start = zero_T, lvar = -max_u_phi, uvar = max_u_phi)
    # Final time
    c, tf = add_var(c, 1; start = 1.0, lvar = zero_T)

    # Expressions (matching JuMP @expressions)
    c, I_the = add_expr(c, ((L - rho[i])^3 + rho[i]^3) * (sin(phi[i]))^2 * third for i = 1:args.nh+1)
    c, I_phi = add_expr(c, ((L - rho[i])^3 + rho[i]^3) * third for i = 1:args.nh+1)

    c, _ = add_obj(c, tf[1] for i = 1:1)

    # Dynamics
    c, _ = add_con(c, -rho[j] + rho[j-1] + half * tf[1] * inv_nh * (rho_dot[j] + rho_dot[j-1]) for j = 2:args.nh+1)
    c, _ = add_con(c, -phi[j] + phi[j-1] + half * tf[1] * inv_nh * (phi_dot[j] + phi_dot[j-1]) for j = 2:args.nh+1)
    c, _ = add_con(c, -the[j] + the[j-1] + half * tf[1] * inv_nh * (the_dot[j] + the_dot[j-1]) for j = 2:args.nh+1)
    c, _ = add_con(c, -rho_dot[j] + rho_dot[j-1] + half * tf[1] * inv_nh * (u_rho[j] + u_rho[j-1]) / L for j = 2:args.nh+1)
    c, _ = add_con(c, -the_dot[j] + the_dot[j-1] + half * tf[1] * inv_nh * (u_the[j] / I_the[j] + u_the[j-1] / I_the[j-1]) for j = 2:args.nh+1)
    c, _ = add_con(c, -phi_dot[j] + phi_dot[j-1] + half * tf[1] * inv_nh * (u_phi[j] / I_phi[j] + u_phi[j-1] / I_phi[j-1]) for j = 2:args.nh+1)

    # Boundary conditions
    c, _ = add_con(c, -rho[1] + rho0 for _ in 1:1)
    c, _ = add_con(c, -the[1] + zero_T for _ in 1:1)
    c, _ = add_con(c, -phi[1] + phi0 for _ in 1:1)
    c, _ = add_con(c, -rho[k] + rho0 for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, -the[k] + two_pi_3 for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, -phi[k] + phi0 for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, -rho_dot[1] + zero_T for _ in 1:1)
    c, _ = add_con(c, -the_dot[1] + zero_T for _ in 1:1)
    c, _ = add_con(c, -phi_dot[1] + zero_T for _ in 1:1)
    c, _ = add_con(c, -rho_dot[k] + zero_T for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, -the_dot[k] + zero_T for k in args.nh+1:args.nh+1)
    c, _ = add_con(c, -phi_dot[k] + zero_T for k in args.nh+1:args.nh+1)
    c
end
