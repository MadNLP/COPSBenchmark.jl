# Rocket Steering Problem
# Trapezoidal formulation
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.steering_model(::JuMPBackend, nh)
    a     = COPSBenchmark.steering_a
    u_min = COPSBenchmark.steering_u_min
    u_max = COPSBenchmark.steering_u_max
    xs    = COPSBenchmark.steering_xs
    xf    = COPSBenchmark.steering_xf

    model = Model()

    @variable(model, u_min <= u[i=1:nh+1] <= u_max, start=0.0)
    @variable(model, x[i=1:nh+1, j=1:4], start=COPSBenchmark.steering_x0(i, j, nh))
    @variable(model, tf, start=1.0)

    @expression(model, h, tf / nh)
    @objective(model, Min, tf)

    @constraint(model, tf >= 0.0)
    @constraints(
        model, begin
            [i=1:nh], x[i+1,1] == x[i,1] + 0.5*h*(x[i,3] + x[i+1,3])
            [i=1:nh], x[i+1,2] == x[i,2] + 0.5*h*(x[i,4] + x[i+1,4])
            [i=1:nh], x[i+1,3] == x[i,3] + 0.5*h*(a*cos(u[i]) + a*cos(u[i+1]))
            [i=1:nh], x[i+1,4] == x[i,4] + 0.5*h*(a*sin(u[i]) + a*sin(u[i+1]))
        end
    )
    @constraint(model, [j=1:4], x[1, j] == xs[j])
    @constraint(model, [j=2:4], x[nh+1, j] == xf[j])
    return model
end
