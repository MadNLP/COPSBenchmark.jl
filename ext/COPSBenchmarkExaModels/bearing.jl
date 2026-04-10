# Journal bearing problem
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.bearing_model(::ExaModelsBackend, nx, ny; T = Float64, backend = nothing, kwargs...)
    b = COPSBenchmark.bearing_b
    e = COPSBenchmark.bearing_e

    hx = 2*pi / (nx+1)
    hy = 2*b / (ny+1)

    v0 = [COPSBenchmark.bearing_v0_ij(i, hx) for i in 1:nx+2, j in 1:ny+2]

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(core, v, 1:nx+2, 1:ny+2; lvar = 0.0, start = v0)

    ExaModels.@add_objective(
        core,
        0.5*(hx*hy/6.0) * (COPSBenchmark.bearing_wq(i, hx) + 2*COPSBenchmark.bearing_wq(i+1, hx))*(((v[i+1,j]-v[i,j])/hx)^2 + ((v[i,j+1]-v[i,j])/hy)^2) for i in 1:nx+1, j in 1:ny+1
    )
    ExaModels.@add_objective(
        core,
        0.5*(hx*hy/6.0) * (2*COPSBenchmark.bearing_wq(i, hx) + 2*COPSBenchmark.bearing_wq(i-1, hx))*(((v[i-1,j]-v[i,j])/hx)^2 + ((v[i,j-1]-v[i,j])/hy)^2) for i in 2:nx+2, j in 2:ny+2
    )
    ExaModels.@add_objective(
        core,
        -hx*hy*e*sin((i-1)*hx)*v[i, j] for i in 1:nx+2, j in 1:ny+2
    )

    ExaModels.@add_constraint(core, c1, v[i, 1] for i in 1:nx+2)
    ExaModels.@add_constraint(core, c2, v[i, ny+2] for i in 1:nx+2)
    ExaModels.@add_constraint(core, c3, v[1, i] for i in 1:ny+2)
    ExaModels.@add_constraint(core, c4, v[nx+2, i] for i in 1:ny+2)

    return ExaModels.ExaModel(core; kwargs...)
end
