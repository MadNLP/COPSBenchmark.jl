# Torsion problem
# Liz Dolan - Summer 2000
# Version 2.0 - October 2000
# COPS 3.1 - March 2004

# Two sizes again.  The spacings and the triangle area keep their symbolic form;
# the distance-to-boundary field, the row list and the bound vectors built from
# it are comprehensions over the grid and travel as data.
@inline function COPSBenchmark.torsion_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    c_val = T(5.0)

    core, nx, ny, d = ExaModels.ExaCore(T; backend = backend, nargs = Val(3))

    hx = one(T) / (nx + 1)
    hy = one(T) / (ny + 1)
    area = T(0.5) * hx * hy

    ExaModels.@add_var(core, v, nx+2, ny+2; start = d.D)

    # Objective = area * ((quadLower + quadUpper)/2 - c*(linLower + linUpper)/3)
    ExaModels.@add_obj(core, area / 2 * (((v[k1+1,k2] - v[k1,k2])/hx)^2 + ((v[k1,k2+1] - v[k1,k2])/hy)^2)
                   for k1 in 1:nx+1, k2 in 1:ny+1)
    ExaModels.@add_obj(core, area / 2 * (((v[k1,k2] - v[k1-1,k2])/hx)^2 + ((v[k1,k2] - v[k1,k2-1])/hy)^2)
                   for k1 in 2:nx+2, k2 in 2:ny+2)
    ExaModels.@add_obj(core, -area * c_val / 3 * (v[k1+1,k2] + v[k1,k2] + v[k1,k2+1])
                   for k1 in 1:nx+1, k2 in 1:ny+1)
    ExaModels.@add_obj(core, -area * c_val / 3 * (v[k1,k2] + v[k1-1,k2] + v[k1,k2-1])
                   for k1 in 2:nx+2, k2 in 2:ny+2)

    # Bound constraints on v: -D <= v <= D (matching JuMP's @constraint formulation)
    ExaModels.@add_con(core, c1, v[k1, k2] for (k1, k2, dd) in d.D_flat; lcon = d.lcon, ucon = d.ucon)

    return core
end

@inline function COPSBenchmark.torsion_args(::ExaModelsBackend, nx, ny; T = Float64)
    hx = T(1.0 / (nx + 1.0))
    hy = T(1.0 / (ny + 1.0))
    D = T[min(min(i, nx-i+1)*hx, min(j, ny-j+1)*hy) for i in 0:nx+1, j in 0:ny+1]
    D_flat = [(k1, k2, D[k1, k2]) for k1 in 1:nx+2, k2 in 1:ny+2]
    lcon = [-d for (_, _, d) in D_flat]
    ucon = [d for (_, _, d) in D_flat]
    return (nx, ny, (; D, D_flat, lcon, ucon))
end

@inline COPSBenchmark.torsion_model(b::ExaModelsBackend, nx, ny; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.torsion_recipe(b; T = T, backend = backend),
        COPSBenchmark.torsion_args(b, nx, ny; T = T)...;
        kwargs...,
    )
