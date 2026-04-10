# Robot Arm Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.robot_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    L         = COPSBenchmark.robot_L
    max_u_rho = COPSBenchmark.robot_max_u_rho
    max_u_the = COPSBenchmark.robot_max_u_the
    max_u_phi = COPSBenchmark.robot_max_u_phi
    rho0      = COPSBenchmark.robot_rho0
    phi0      = COPSBenchmark.robot_phi0

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(core, rho, nh+1; start=rho0, lvar = 0, uvar = L)
    ExaModels.@add_variable(core, the, nh+1; start=[2*pi/3*(k/nh)^2 for k=1:nh+1], lvar = -pi, uvar = pi)
    ExaModels.@add_variable(core, phi, nh+1; start=phi0, lvar = 0, uvar = pi)
    ExaModels.@add_variable(core, rho_dot, nh+1; start=0.0)
    ExaModels.@add_variable(core, the_dot, nh+1; start=[4*pi/3*(k/nh) for k=1:nh+1])
    ExaModels.@add_variable(core, phi_dot, nh+1; start=0.0)
    ExaModels.@add_variable(core, u_rho, nh+1; start=0.0, lvar = -max_u_rho, uvar = max_u_rho)
    ExaModels.@add_variable(core, u_the, nh+1; start=0.0, lvar = -max_u_the, uvar = max_u_the)
    ExaModels.@add_variable(core, u_phi, nh+1; start=0.0, lvar = -max_u_phi, uvar = max_u_phi)
    ExaModels.@add_variable(core, tf, 1; start = 1.0, lvar = 0.0)

    ExaModels.@add_expression(core, I_the, ((L-rho[i])^3+rho[i]^3)*(sin(phi[i]))^2/3.0 for i=1:nh+1)
    ExaModels.@add_expression(core, I_phi, ((L-rho[i])^3+rho[i]^3)/3.0 for i=1:nh+1)

    ExaModels.@add_objective(core, tf[1] for i=1:1)

    ExaModels.@add_constraint(core, c1, - rho[j] + rho[j-1] + 0.5 * tf[1]/nh * (rho_dot[j] + rho_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_constraint(core, c2, - phi[j] + phi[j-1] + 0.5 * tf[1]/nh * (phi_dot[j] + phi_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_constraint(core, c3, - the[j] + the[j-1] + 0.5 * tf[1]/nh * (the_dot[j] + the_dot[j-1]) for j=2:nh+1)
    ExaModels.@add_constraint(core, c4, - rho_dot[j] + rho_dot[j-1] + 0.5 * tf[1]/nh * (u_rho[j] + u_rho[j-1]) / L for j=2:nh+1)
    ExaModels.@add_constraint(core, c5, - the_dot[j] + the_dot[j-1] + 0.5 * tf[1]/nh * (u_the[j] / I_the[j] + u_the[j-1] / I_the[j-1]) for j=2:nh+1)
    ExaModels.@add_constraint(core, c6, - phi_dot[j] + phi_dot[j-1] + 0.5 * tf[1]/nh * (u_phi[j] / I_phi[j] + u_phi[j-1] / I_phi[j-1]) for j=2:nh+1)

    ExaModels.@add_constraint(core, c7,  - rho[1] + 4.5)
    ExaModels.@add_constraint(core, c8,  - the[1] + 0.0)
    ExaModels.@add_constraint(core, c9,  - phi[1] + pi / 4.0)
    ExaModels.@add_constraint(core, c10, - rho[nh+1] + 4.5)
    ExaModels.@add_constraint(core, c11, - the[nh+1] + 2.0 * pi / 3)
    ExaModels.@add_constraint(core, c12, - phi[nh+1] + pi / 4.0)
    ExaModels.@add_constraint(core, c13, - rho_dot[1] + 0.0)
    ExaModels.@add_constraint(core, c14, - the_dot[1] + 0.0)
    ExaModels.@add_constraint(core, c15, - phi_dot[1] + 0.0)
    ExaModels.@add_constraint(core, c16, - rho_dot[nh+1] + 0.0)
    ExaModels.@add_constraint(core, c17, - the_dot[nh+1] + 0.0)
    ExaModels.@add_constraint(core, c18, - phi_dot[nh+1] + 0.0)

    return ExaModels.ExaModel(core; kwargs...)
end
