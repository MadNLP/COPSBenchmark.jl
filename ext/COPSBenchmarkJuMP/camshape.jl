# Cam Shape Problem
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.camshape_model(::JuMPBackend, n)
    R_v   = COPSBenchmark.camshape_R_v
    R_max = COPSBenchmark.camshape_R_max
    R_min = COPSBenchmark.camshape_R_min
    alpha = COPSBenchmark.camshape_alpha

    d_theta = COPSBenchmark.camshape_d_theta(n)

    model = Model()

    @variable(model, R_min <= r[1:n] <= R_max, start=(R_min+R_max)/2.0)

    @objective(model, Max, ((pi*R_v)/n) * sum(r[i] for i in 1:n))

    @constraint(
        model,
        [i=2:n-1],
        COPSBenchmark.camshape_convexity(r[i-1], r[i], r[i+1], d_theta) <= 0.0
    )
    @constraints(
        model, begin
            COPSBenchmark.camshape_convexity(R_min, r[1], r[2], d_theta) <= 0.0
            - R_min^2 - R_min*r[1] + 2*R_min*r[1]*cos(d_theta) <= 0.0
            COPSBenchmark.camshape_convexity(r[n-1], r[n], R_max, d_theta) <= 0.0
            - 2*R_max*r[n] + 2*r[n]^2*cos(d_theta) <= 0.0
        end
    )
    @constraint(
        model,
        [i=1:n-1],
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
