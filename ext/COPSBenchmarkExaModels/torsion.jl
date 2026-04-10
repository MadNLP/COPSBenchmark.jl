# Torsion problem
# Liz Dolan - Summer 2000
# Version 2.0 - October 2000
# COPS 3.1 - March 2004

@inline function COPSBenchmark.torsion_model(::ExaModelsBackend, nx, ny; T = Float64, backend = nothing, kwargs...)
    c_val = T(5.0)
    hx = T(1.0 / (nx + 1.0))
    hy = T(1.0 / (ny + 1.0))
    area = T(0.5) * hx * hy

    # Distance to boundary: D[k1,k2] for k1 in 1:nx+2, k2 in 1:ny+2
    D = T[min(min(i, nx-i+1)*hx, min(j, ny-j+1)*hy) for i in 0:nx+1, j in 0:ny+1]

    core = ExaModels.ExaCore(T; backend = backend)
    ExaModels.@add_variable(core, v, nx+2, ny+2; start = D)

    # Objective = area * ((quadLower + quadUpper)/2 - c*(linLower + linUpper)/3)
    # quadLower: i in 0:nx, j in 0:ny → k1 in 1:nx+1, k2 in 1:ny+1
    ExaModels.@add_objective(core, area / 2 * (((v[k1+1,k2] - v[k1,k2])/hx)^2 + ((v[k1,k2+1] - v[k1,k2])/hy)^2)
                   for k1 in 1:nx+1, k2 in 1:ny+1)
    # quadUpper: i in 1:nx+1, j in 1:ny+1 → k1 in 2:nx+2, k2 in 2:ny+2
    ExaModels.@add_objective(core, area / 2 * (((v[k1,k2] - v[k1-1,k2])/hx)^2 + ((v[k1,k2] - v[k1,k2-1])/hy)^2)
                   for k1 in 2:nx+2, k2 in 2:ny+2)
    # linLower: i in 0:nx, j in 0:ny → k1 in 1:nx+1, k2 in 1:ny+1
    ExaModels.@add_objective(core, -area * c_val / 3 * (v[k1+1,k2] + v[k1,k2] + v[k1,k2+1])
                   for k1 in 1:nx+1, k2 in 1:ny+1)
    # linUpper: i in 1:nx+1, j in 1:ny+1 → k1 in 2:nx+2, k2 in 2:ny+2
    ExaModels.@add_objective(core, -area * c_val / 3 * (v[k1,k2] + v[k1-1,k2] + v[k1,k2-1])
                   for k1 in 2:nx+2, k2 in 2:ny+2)

    # Bound constraints on v: -D <= v <= D (matching JuMP's @constraint formulation)
    D_flat = [(k1, k2, D[k1, k2]) for k1 in 1:nx+2, k2 in 1:ny+2]
    ExaModels.@add_constraint(core, c1, v[k1, k2] for (k1, k2, d) in D_flat; lcon = [-d for (_, _, d) in D_flat], ucon = [d for (_, _, d) in D_flat])

    return ExaModels.ExaModel(core; kwargs...)
end
