# Robot Arm Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function robot_model(nh)
    # total length of arm
    L = 5.0

    # Upper bounds on the controls
    max_u_rho = 1.0
    max_u_the = 1.0
    max_u_phi = 1.0
    # Initial positions of the length and the angles for the robot arm
    rho0 = 4.5
    phi0 = pi /4

    model = Model()

    @variable(model, 0 <= rho[k=1:nh+1] <= L, start=rho0)
    @variable(model, -pi <= the[k=1:nh+1] <= pi, start=2*pi/3*(k/nh)^2)
    @variable(model, 0 <= phi[k=1:nh+1] <= pi, start=phi0)
    # Derivatives
    @variable(model, rho_dot[k=1:nh+1], start=0.0)
    @variable(model, the_dot[k=1:nh+1], start=4*pi/3*(k/nh))
    @variable(model, phi_dot[k=1:nh+1], start=0.0)
    # Control
    @variable(model, -max_u_rho <= u_rho[1:nh+1] <= max_u_rho, start=0.0)
    @variable(model, -max_u_the <= u_the[1:nh+1] <= max_u_the, start=0.0)
    @variable(model, -max_u_phi <= u_phi[1:nh+1] <= max_u_phi, start=0.0)
    # Steps and final time
    @variable(model, step >= 0.0)
    # The moments of inertia
    @variable(model, I_the[1:nh+1], start=((L-rho0)^3+rho0^3)*(sin(phi0))^2/3.0)
    @variable(model, I_phi[1:nh+1], start=((L-rho0)^3+rho0^3)/3.0)

    @objective(model, Min, step * nh)

    # Physical equations
    @constraint(
        model,
        [i=1:nh+1],
        I_the[i] == ((L-rho[i])^3+rho[i]^3)*(sin(phi[i]))^2/3.0
    )
    @constraint(
        model,
        [i=1:nh+1],
        I_phi[i] == ((L-rho[i])^3+rho[i]^3)/3.0
    )
    # Dynamics
    @constraint(
        model,
        [j=2:nh+1],
        rho[j] == rho[j-1] + 0.5 * step * (rho_dot[j] + rho_dot[j-1]),
    )
    @constraint(
        model,
        [j=2:nh+1],
        phi[j] == phi[j-1] + 0.5 * step * (phi_dot[j] + phi_dot[j-1]),
    )
    @constraint(
        model,
        [j=2:nh+1],
        the[j] == the[j-1] + 0.5 * step * (the_dot[j] + the_dot[j-1]),
    )

    @constraint(
        model,
        [j=2:nh+1],
        rho_dot[j] == rho_dot[j-1] + 0.5 * step * (u_rho[j] + u_rho[j-1]) / L
    )
    @constraint(
        model,
        [j=2:nh+1],
        the_dot[j] == the_dot[j-1] + 0.5 * step * (u_the[j] / I_the[j] + u_the[j-1] / I_the[j-1])
    )
    @constraint(
        model,
        [j=2:nh+1],
        phi_dot[j] == phi_dot[j-1] + 0.5 * step * (u_phi[j] / I_phi[j] + u_phi[j-1] / I_phi[j-1])
    )
    # Boundary condition
    @constraints(
        model, begin
            rho[1] == 4.5
            the[1] == 0.0
            phi[1] == pi / 4.0
            rho[nh+1] == 4.5
            the[nh+1] == 2.0 * pi / 3
            phi[nh+1] == pi / 4.0
            rho_dot[1] == 0.0
            the_dot[1] == 0.0
            phi_dot[1] == 0.0
            rho_dot[nh+1] == 0.0
            the_dot[nh+1] == 0.0
            phi_dot[nh+1] == 0.0
        end
    )

    return model
end

