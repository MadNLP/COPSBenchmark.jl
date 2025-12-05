# Robot Arm Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.robot_model(nh, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)

    # total length of arm
    L = 5.0

    # Upper bounds on the controls
    max_u_rho = 1.0
    max_u_the = 1.0
    max_u_phi = 1.0
    # Initial positions of the length and the angles for the robot arm
    rho0 = 4.5
    phi0 = pi /4

    core = ExaModels.ExaCore(T; backend= backend)

    rho = ExaModels.variable(core,nh+1; start=rho0, lvar =  0 , uvar = L)
    the = ExaModels.variable(core, nh+1; start=[2*pi/3*(k/nh)^2 for k=1:nh+1], lvar =  -pi , uvar = pi)
    phi = ExaModels.variable(core, nh+1; start=phi0, lvar = 0, uvar = pi)
    # Derivatives
    rho_dot = ExaModels.variable(core, nh+1; start=0.0)
    the_dot = ExaModels.variable(core, nh+1; start=[4*pi/3*(k/nh) for k=1:nh+1])
    phi_dot = ExaModels.variable(core, nh+1; start=0.0)
    # Control
    u_rho = ExaModels.variable(core, nh+1; start=0.0, lvar = -max_u_rho, uvar = max_u_rho)
    u_the = ExaModels.variable(core, nh+1; start=0.0, lvar = -max_u_the, uvar = max_u_the)
    u_phi = ExaModels.variable(core, nh+1; start=0.0, lvar = -max_u_phi, uvar = max_u_phi)
    # Steps and final time
    step = ExaModels.variable(core, lvar = 0.0)
    # The moments of inertia
    I_the = ExaModels.variable(core, nh+1; start=((L-rho0)^3+rho0^3)*(sin(phi0))^2/3.0)
    I_phi = ExaModels.variable(core, nh+1; start=((L-rho0)^3+rho0^3)/3.0)

    ExaModels.objective(core, step * nh for i=1:1)

    # Physical equations
    ExaModels.constraint(
        core,
        - I_the[i] + ((L-rho[i])^3+rho[i]^3)*(sin(phi[i]))^2/3.0 for i=1:nh+1
    )
    ExaModels.constraint(
        core,
        - I_phi[i]+ ((L-rho[i])^3+rho[i]^3)/3.0 for i=1:nh+1
    )
    # Dynamics
    ExaModels.constraint(
        core,
        - rho[j] + rho[j-1] + 0.5 * step * (rho_dot[j] + rho_dot[j-1]) for j=2:nh+1
            )
    ExaModels.constraint(
        core,
        - phi[j] + phi[j-1] + 0.5 * step * (phi_dot[j] + phi_dot[j-1]) for j=2:nh+1
            )
    ExaModels.constraint(
        core,
        - the[j] + the[j-1] + 0.5 * step * (the_dot[j] + the_dot[j-1]) for j=2:nh+1
    )

    ExaModels.constraint(
        core,
        - rho_dot[j] + rho_dot[j-1] + 0.5 * step * (u_rho[j] + u_rho[j-1]) / L for j=2:nh+1
    )
    ExaModels.constraint(
        core,
        - the_dot[j] + the_dot[j-1] + 0.5 * step * (u_the[j] / I_the[j] + u_the[j-1] / I_the[j-1]) for j=2:nh+1
    )
    c = ExaModels.constraint(
        core,
        - phi_dot[j] + phi_dot[j-1] + 0.5 * step * (u_phi[j] / I_phi[j] + u_phi[j-1] / I_phi[j-1]) for j=2:nh+1
            )
    ExaModels.constraint(
        core, - rho[1] + 4.5
    )
    ExaModels.constraint(
        core, - the[1] + 0.0
    )
    ExaModels.constraint(
        core, - phi[1] + pi / 4.0
    )
    ExaModels.constraint(
        core, - rho[nh+1] + 4.5
    )
    ExaModels.constraint(
        core, - the[nh+1] + 2.0 * pi / 3
    )
    ExaModels.constraint(
        core, - phi[nh+1] + pi / 4.0
    )
    ExaModels.constraint(
        core, - rho_dot[1] + 0.0
    )
    ExaModels.constraint(
        core, - the_dot[1] + 0.0
    )
    ExaModels.constraint(
        core, - phi_dot[1] + 0.0
    )
    ExaModels.constraint(
        core, - rho_dot[nh+1] + 0.0
    )
    ExaModels.constraint(
        core, - the_dot[nh+1] + 0.0
    )
    ExaModels.constraint(
        core, - phi_dot[nh+1] + 0.0
    )

    return ExaModels.ExaModel(core; kwargs...)
end

