# Journal bearing problem
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function bearing_model(nx, ny)
    b = 10              # grid is (0,2*pi)x(0,2*b)
    e = 0.1             # eccentricity

    hx = 2*pi / (nx+1)  # grid spacing
    hy = 2*b / (ny+1)   # grid spacing
    area = 0.5*hx*hy    # area of triangle

    wq = [(1.0 + e*cos(i*hx))^3 for i in 0:nx+1]

    model = Model()

    @variable(model, v[i=1:nx+2, j=1:ny+2] >= 0.0, start=max(sin((i-1)*hx), 0.0))

    @expression(
        model,
        ex1,
        sum((wq[i] + 2*wq[i+1])*(((v[i+1,j]-v[i,j])/hx)^2 + ((v[i,j+1]-v[i,j])/hy)^2) for i in 1:nx+1, j in 1:ny+1)
    )
    @expression(
        model,
        ex2,
        sum((2*wq[i] + 2*wq[i-1])*(((v[i-1,j]-v[i,j])/hx)^2 + ((v[i,j-1]-v[i,j])/hy)^2) for i in 2:nx+2, j in 2:ny+2)
    )
    @expression(
        model,
        ex3,
        sum(e*sin((i-1)*hx)*v[i, j] for i in 1:nx+2, j in 1:ny+2),
    )

    @objective(model, Min, 0.5*(hx*hy/6.0) * (ex1 + ex2) - hx*hy*ex3)
    # Boundary condition
    @constraint(model, [i=1:nx+2], v[i, 1] == 0.0)
    @constraint(model, [i=1:nx+2], v[i, ny+2] == 0.0)
    @constraint(model, [j=1:nx+2], v[1, j] == 0.0)
    @constraint(model, [j=1:nx+2], v[nx+2, j] == 0.0)

    return model
end

