# Journal bearing problem
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.bearing_model(nx, ny, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)
    b = 10              # grid is (0,2*pi)x(0,2*b)
    e = 0.1             # eccentricity

    hx = 2*pi / (nx+1)  # grid spacing
    hy = 2*b / (ny+1)   # grid spacing
    area = 0.5*hx*hy    # area of triangle

    wq(i) = (1.0 + e*cos((i-1)*hx))^3
    v0 = [max(sin((i-1)*hx), 0.0) for i in 1:nx+2, j in 1:ny+2]

    core = ExaModels.ExaCore(T; backend=backend)

    v = ExaModels.variable(core, 1:nx+2, 1:ny+2; lvar = 0.0, start=v0)

    ExaModels.objective(
        core,
        0.5*(hx*hy/6.0) * (wq(i) + 2*wq(i+1))*(((v[i+1,j]-v[i,j])/hx)^2 + ((v[i,j+1]-v[i,j])/hy)^2) for i in 1:nx+1, j in 1:ny+1
    )
    ExaModels.objective(
        core,
        0.5*(hx*hy/6.0) * (2*wq(i) + 2*wq(i-1))*(((v[i-1,j]-v[i,j])/hx)^2 + ((v[i,j-1]-v[i,j])/hy)^2) for i in 2:nx+2, j in 2:ny+2
    )
    ExaModels.objective(
        core,
        -hx*hy*e*sin((i-1)*hx)*v[i, j] for i in 1:nx+2, j in 1:ny+2
    )

    ExaModels.constraint(core, v[i, 1] for i in 1:nx+2)
    ExaModels.constraint(core, v[i, ny+2] for i in 1:nx+2)
    ExaModels.constraint(core, v[1, i] for i in 1:ny+2)
    ExaModels.constraint(core, v[nx+2, i] for i in 1:ny+2)

    return ExaModels.ExaModel(core; kwargs...)
end



