# Flow in a Channel Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

# Every table here is a function of the interval length h = tf/nh, so all of it
# is data rather than structure.  It travels as a single named tuple rather than
# as ten separate arguments: a placeholder supports field access, so the recipe
# reads `d.con1_itr` exactly as the original read `con1_itr`.
@inline function COPSBenchmark.channel_recipe(
    ::ExaModelsBackend; T = Float64, backend = nothing,
)
    bc = T[0.0 1.0; 0.0 0.0]

    core, nh = ExaModels.ExaCore(T; backend = backend, nargs = Val(1))
    d = ExaModels.ArgNode1(COPSBenchmark.channel_data, nh)

    ExaModels.@add_var(core, v, nh, 4; start = d.v0)
    ExaModels.@add_var(core, w, nh, 4; start = 0.0)
    ExaModels.@add_var(core, uc, nh, 4, 4; start = d.uc0)
    ExaModels.@add_var(core, Duc, nh, 4, 4; start = 0.0)

    # Constant objective
    ExaModels.@add_obj(core, one(T) for _i in 1:1)

    # con1: uc[i,j,s] - v[i,s] - h * sum(w[i,k] * rho[j]^k / k!)
    ExaModels.@add_con(core, c1,
        uc[i, j, s] - v[i, s] - (a1*w[i,1] + a2*w[i,2] + a3*w[i,3] + a4*w[i,4])
        for (i, j, s, a1, a2, a3, a4) in d.con1_itr
    )

    # con2: Duc[i,j,s] - sum(v[i,k]*Bv[k]) - sum(w[i,k]*Bw[k])
    ExaModels.@add_con(core, c2,
        Duc[i, j, s] - (bv1*v[i,1] + bv2*v[i,2] + bv3*v[i,3] + bv4*v[i,4])
                     - (bw1*w[i,1] + bw2*w[i,2] + bw3*w[i,3] + bw4*w[i,4])
        for (i, j, s, bv1, bv2, bv3, bv4, bw1, bw2, bw3, bw4) in d.con2_itr
    )

    # Boundary conditions.  bc3/bc4 index the last interval, so the index comes
    # from a one-element set and the coefficients ride along with it.
    ExaModels.@add_con(core, bc1, v[1, 1] - bc[1, 1])
    ExaModels.@add_con(core, bc2, v[1, 2] - bc[2, 1])

    ExaModels.@add_con(core, bc3,
        p1*v[k,1] + p2*v[k,2] + p3*v[k,3] + p4*v[k,4] +
        q1*w[k,1] + q2*w[k,2] + q3*w[k,3] + q4*w[k,4] - bc[1, 2]
        for (k, p1, p2, p3, p4, q1, q2, q3, q4) in d.bc3_itr
    )

    ExaModels.@add_con(core, bc4,
        p1*v[k,2] + p2*v[k,3] + p3*v[k,4] +
        q1*w[k,1] + q2*w[k,2] + q3*w[k,3] + q4*w[k,4] - bc[2, 2]
        for (k, p1, p2, p3, q1, q2, q3, q4) in d.bc4_itr
    )

    # Continuity
    ExaModels.@add_con(core, c3,
        cv1*v[i,1] + cv2*v[i,2] + cv3*v[i,3] + cv4*v[i,4] +
        cw1*w[i,1] + cw2*w[i,2] + cw3*w[i,3] + cw4*w[i,4] - v[i+1, s]
        for (i, s, cv1, cv2, cv3, cv4, cw1, cw2, cw3, cw4) in d.cont_itr
    )

    # Collocation physics
    ExaModels.@add_con(core, c4,
        e1*w[i,1] + e2*w[i,2] + e3*w[i,3] + e4*w[i,4] -
        T(10.0) * (Duc[i, j, 2] * Duc[i, j, 3] - Duc[i, j, 1] * Duc[i, j, 4])
        for (i, j, e1, e2, e3, e4) in d.coll_itr
    )

    return core
end

@inline COPSBenchmark.channel_args(::ExaModelsBackend, nh; T = Float64) = (nh,)

@inline COPSBenchmark.channel_model(b::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...) =
    ExaModels.ExaModel(
        COPSBenchmark.channel_recipe(b; T = T, backend = backend),
        COPSBenchmark.channel_args(b, nh; T = T)...;
        kwargs...,
    )
