# Cam Shape Problem
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.camshape_model(n, ::JuMPBackend)
    R_v = 1.0         # design parameter related to the valve shape
    R_max = 2.0       # maximum allowed radius of the cam
    R_min = 1.0       # minimum allowed radius of the cam
    alpha = 1.5       # curvature limit parameter

    d_theta = 2*pi/(5*(n+1))   # angle between discretization points

    model = Model()
    # radius of the cam at discretization points
    @variable(model, R_min <= r[1:n] <= R_max, start=(R_min+R_max)/2.0)

    @objective(model, Max, ((pi*R_v)/n) * sum(r[i] for i in 1:n))

    # Convexity
    @constraint(
        model,
        [i=2:n-1],
        - r[i-1]*r[i] - r[i]*r[i+1] + 2*r[i-1]*r[i+1]*cos(d_theta) <= 0.0
    )
    @constraints(
        model, begin
            - R_min*r[1] - r[1]*r[2] + 2*R_min*r[2]*cos(d_theta) <= 0.0
            - R_min^2 - R_min*r[1] + 2*R_min*r[1]*cos(d_theta) <= 0.0
            - r[n-1]*r[n] - r[n]*R_max + 2*r[n-1]*R_max*cos(d_theta) <= 0.0
            - 2*R_max*r[n] + 2*r[n]^2*cos(d_theta) <= 0.0
        end
    )
    # Curvature
    @constraint(
        model,
        [i=1:n-1] ,
        -alpha*d_theta <= (r[i+1] - r[i]) <= alpha*d_theta,
    )
    @constraints(
        model, begin
            -alpha*d_theta <= (r[1] - R_min) <= alpha*d_theta
            -alpha*d_theta <= (R_max - r[n]) <= alpha*d_theta
        end
    )

    return model
end

