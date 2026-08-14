# Minimal surface with obstacle problem

#  Find the surface with minimal area, given boundary conditions,
#  and above an obstacle.

#  This is problem 17=the COPS (Version 3) collection of
#  E. Dolan and J. More'
#  see "Benchmarking Optimization Software with COPS"
#  Argonne National Labs Technical Report ANL/MCS-246 (2004)
#  classification OBR2-AN-V-V

# The spacings and the triangle area stay symbolic.  The surface initialisation
# is a comprehension over the mesh, and the interior block's index ranges are
# `Int(floor(0.25/hx))`-style constructor calls, which a placeholder cannot take
# -- both are data.
@inline function COPSBenchmark.minsurf_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    c, nx, ny = ExaModels.ExaCore(T; backend = backend, nargs = Val(2))
    d = ExaModels.ArgNode2(COPSBenchmark.minsurf_data, nx, ny)

    hx = 1 / (nx + 1)
    hy = 1 / (ny + 1)
    area = 1 // 2 * hx * hy

    ExaModels.@add_var(c, v, nx+2, ny+2; start = d.v0)

    ExaModels.@add_obj(c, area * (1 + ((v[i + 1, j] - v[i, j]) / hx)^2 + ((v[i, j + 1] - v[i, j]) / hy)^2)^(1 / 2) for
                        i = 1:(nx + 1), j = 1:(ny + 1))
    ExaModels.@add_obj(c, area * (1 + ((v[i - 1, j] - v[i, j]) / hx)^2 + ((v[i, j - 1] - v[i, j]) / hy)^2)^(1 / 2) for
                         i = 2:(nx + 2), j = 2:(ny + 2))

    ExaModels.@add_con(
        c,
        c1,
        v[1, j + 1] for j in 0:ny+1
    )

    # The fixed edge indices depend on the size, so they come from one-element
    # sets; the enclosing set has a single value, so the row order is unchanged.
    ExaModels.@add_con(
        c,
        c2,
        v[k, j + 1] for k in (nx+2):(nx+2), j in 0:ny+1
    )

    ExaModels.@add_con(
        c,
        c3,
        v[i + 1, 1] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1
    )

    ExaModels.@add_con(
        c,
        c4,
        v[i + 1, k] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1, k in (ny+2):(ny+2)
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
        v[i + 1, j + 1] for i in d.c6_i, j in d.c6_j;
        lcon = 1,
        ucon = Inf
    )

    return c
end

@inline COPSBenchmark.minsurf_args(::ExaModelsBackend, nx::Int, ny::Int) = (nx, ny)

@inline COPSBenchmark.minsurf_model(b::ExaModelsBackend, nx::Int, ny::Int; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.minsurf_recipe(b; T = T, backend = backend),
        COPSBenchmark.minsurf_args(b, nx, ny)...;
        kwargs...,
    )
