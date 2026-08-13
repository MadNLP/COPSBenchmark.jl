# Journal bearing problem
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# Two sizes, so two placeholders.  The grid spacings depend on them and keep
# their symbolic form; the starting field and the per-column weights are
# comprehensions over the grid and travel as data.
@inline function COPSBenchmark.bearing_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    b = 10              # grid is (0,2*pi)x(0,2*b)
    e = 0.1             # eccentricity

    core, nx, ny, d = ExaModels.ExaCore(T; backend = backend, nargs = Val(3))

    hx = 2*pi / (nx+1)  # grid spacing
    hy = 2*b / (ny+1)   # grid spacing

    ExaModels.@add_var(core, v, 1:nx+2, 1:ny+2; lvar = 0.0, start = d.v0)

    ExaModels.@add_obj(
        core,
        0.5*(hx*hy/6.0) * (w1 + 2*w2)*(((v[i+1,j]-v[i,j])/hx)^2 + ((v[i,j+1]-v[i,j])/hy)^2) for (i, j, w1, w2) in d.lower
    )
    ExaModels.@add_obj(
        core,
        0.5*(hx*hy/6.0) * (2*w1 + 2*w2)*(((v[i-1,j]-v[i,j])/hx)^2 + ((v[i,j-1]-v[i,j])/hy)^2) for (i, j, w1, w2) in d.upper
    )
    ExaModels.@add_obj(
        core,
        -hx*hy*e*s*v[i, j] for (i, j, s) in d.lin
    )

    ExaModels.@add_con(core, c1, v[i, 1] for i in 1:nx+2)
    ExaModels.@add_con(core, c2, v[i, k] for i in 1:nx+2, k in (ny+2):(ny+2))
    ExaModels.@add_con(core, c3, v[1, i] for i in 1:ny+2)
    ExaModels.@add_con(core, c4, v[k, i] for k in (nx+2):(nx+2), i in 1:ny+2)

    return core
end

@inline function COPSBenchmark.bearing_args(::ExaModelsBackend, nx, ny; T = Float64)
    b = 10
    e = 0.1
    hx = 2*pi / (nx+1)
    hy = 2*b / (ny+1)
    wq(i) = (1.0 + e*cos((i-1)*hx))^3
    v0 = [max(sin((i-1)*hx), 0.0) for i in 1:nx+2, j in 1:ny+2]
    # `wq` and the eccentricity term are per-row coefficients; they ride with
    # the indices they belong to rather than being recomputed symbolically.
    lower = [(i, j, wq(i), wq(i+1)) for i in 1:nx+1, j in 1:ny+1]
    upper = [(i, j, wq(i), wq(i-1)) for i in 2:nx+2, j in 2:ny+2]
    lin   = [(i, j, sin((i-1)*hx)) for i in 1:nx+2, j in 1:ny+2]
    return (nx, ny, (; v0, lower, upper, lin))
end

@inline COPSBenchmark.bearing_model(b::ExaModelsBackend, nx, ny; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.bearing_recipe(b; T = T, backend = backend),
        COPSBenchmark.bearing_args(b, nx, ny; T = T)...;
        kwargs...,
    )
