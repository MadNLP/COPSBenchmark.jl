# Minimal surface with obstacle problem

#  Find the surface with minimal area, given boundary conditions,
#  and above an obstacle.

#  This is problem 17=the COPS (Version 3) collection of
#  E. Dolan and J. More'
#  see "Benchmarking Optimization Software with COPS"
#  Argonne National Labs Technical Report ANL/MCS-246 (2004)
#  classification OBR2-AN-V-V

@inline function COPSBenchmark.minsurf_model(::ExaModelsBackend, nx::Int, ny::Int; T = Float64, backend = nothing, kwargs...)
    x_mesh = LinRange(0, 1, nx + 2) # coordinates of the mesh points x

    v0 = zeros(nx + 2, ny + 2) # Surface matrix initialization
    for i = 1:(nx + 2), j = 1:(ny + 2)
        v0[i, j] = 1 - (2 * x_mesh[i] - 1)^2
    end

    hx = 1 / (nx + 1)
    hy = 1 / (ny + 1)
    area = 1 // 2 * hx * hy

    c = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))
    ExaModels.@add_var(c, v, nx+2, ny+2; start = v0)

    ExaModels.@add_obj(c, area * (1 + ((v[i + 1, j] - v[i, j]) / hx)^2 + ((v[i, j + 1] - v[i, j]) / hy)^2)^(1 / 2) for
                        i = 1:(nx + 1), j = 1:(ny + 1))
    ExaModels.@add_obj(c, area * (1 + ((v[i - 1, j] - v[i, j]) / hx)^2 + ((v[i, j - 1] - v[i, j]) / hy)^2)^(1 / 2) for
                         i = 2:(nx + 2), j = 2:(ny + 2))

    ExaModels.@add_con(
        c,
        c1,
        v[1, j + 1] for j in 0:ny+1
    )

    ExaModels.@add_con(
        c,
        c2,
        v[nx + 2, j + 1] for j in 0:ny+1
    )

    ExaModels.@add_con(
        c,
        c3,
        v[i + 1, 1] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1
    )

    ExaModels.@add_con(
        c,
        c4,
        v[i + 1, ny+2] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1
    )

    ExaModels.@add_con(
        c,
        c5,
        v[i + 1, j + 1] for j in 0:ny+1, i in 0:nx+1;
        lcon = 0,
        ucon = Inf
    )

    ExaModels.@add_con(
        c,
        c6,
        v[i + 1, j + 1] for i in Int(floor(0.25 / hx)):Int(ceil(0.75 / hx)), j in Int(floor(0.25 / hy)):Int(ceil(0.75 / hy));
        lcon = 1,
        ucon = Inf
    )

    return ExaModels.ExaModel(c; kwargs...)
end



function COPSBenchmark.minsurf_core()
    args = ExaModels.ArgTracer()
    v0 = Deferred(a -> begin
        x_mesh = LinRange(0, 1, a.nx + 2)
        v0 = zeros(a.nx + 2, a.ny + 2)
        for i = 1:(a.nx + 2), j = 1:(a.ny + 2)
            v0[i, j] = 1 - (2 * x_mesh[i] - 1)^2
        end
        v0
    end)

    hx = 1 / (args.nx + 1)
    hy = 1 / (args.ny + 1)
    area = 1 // 2 * hx * hy

    c = ExaCore(concrete = Val(true))
    c, v = add_var(c, args.nx+2, args.ny+2; start = v0)

    c, _ = add_obj(c, area * (1 + ((v[i + 1, j] - v[i, j]) / hx)^2 + ((v[i, j + 1] - v[i, j]) / hy)^2)^(1 / 2) for
        i = 1:args.nx+1, j = 1:args.ny+1)
    c, _ = add_obj(c, area * (1 + ((v[i - 1, j] - v[i, j]) / hx)^2 + ((v[i, j - 1] - v[i, j]) / hy)^2)^(1 / 2) for
        i = 2:args.nx+2, j = 2:args.ny+2)

    c, _ = add_con(c, v[1, j + 1] for j in 0:args.ny+1)
    c, _ = add_con(c, v[args.nx + 2, j + 1] for j in 0:args.ny+1)
    c, _ = add_con(c, v[i + 1, 1] - 1 + (2 * i * hx - 1)^2 for i in 0:args.nx+1)
    c, _ = add_con(c, v[i + 1, args.ny+2] - 1 + (2 * i * hx - 1)^2 for i in 0:args.nx+1)
    c, _ = add_con(c, (v[i + 1, j + 1] for j in 0:args.ny+1, i in 0:args.nx+1);
        lcon = 0, ucon = Inf)
    c6_itr = Deferred(a -> begin
        hx = 1 / (a.nx + 1)
        hy = 1 / (a.ny + 1)
        [(i = i, j = j) for i in Int(floor(0.25 / hx)):Int(ceil(0.75 / hx)),
                            j in Int(floor(0.25 / hy)):Int(ceil(0.75 / hy))]
    end)
    c, _ = add_con(c, (v[p.i + 1, p.j + 1] for p in c6_itr); lcon = 1, ucon = Inf)
    c
end
