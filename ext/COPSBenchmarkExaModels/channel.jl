# Flow in a Channel Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

@inline function COPSBenchmark.channel_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
    nc = 4
    nd = 4
    R = T(10.0)
    tf = T(1.0)
    h = tf / nh

    bc = T[0.0 1.0; 0.0 0.0]
    rho = T[0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    t = T[(i-1)*h for i in 1:nh+1]

    # Initial value
    v0 = zeros(T, nh, nd)
    for i in 1:nh
        v0[i, 1] = t[i]^2*(3 - 2*t[i])
        v0[i, 2] = 6*t[i]*(1 - t[i])
        v0[i, 3] = 6*(1 - 2*t[i])
        v0[i, 4] = -12
    end
    uc0 = T[v0[i, s] for i in 1:nh, j in 1:nc, s in 1:nd]

    # fac[k+1] = k!
    fac = T[factorial(k) for k in 0:nc+nd]

    # Precompute all coefficients for the collocation constraints
    # con1: uc[i,j,s] = v[i,s] + h * sum_{k=1}^{nc}(w[i,k] * rho[j]^k / k!)
    con1_itr = [
        (i, j, s, h*rho[j]^1/fac[2], h*rho[j]^2/fac[3], h*rho[j]^3/fac[4], h*rho[j]^4/fac[5])
        for i in 1:nh, j in 1:nc, s in 1:nd
    ]

    # con2: Duc[i,j,s] coefficients depend on (j,s)
    con2_itr = [
        (i, j, s,
         # v coefficients: (rho[j]*h)^(k-s) / (k-s)! for k = s:nd (padded to length 4)
         (s <= 1 ? (rho[j]*h)^(1-s)/fac[1-s+1] : zero(T)),
         (s <= 2 ? (rho[j]*h)^(2-s)/fac[2-s+1] : zero(T)),
         (s <= 3 ? (rho[j]*h)^(3-s)/fac[3-s+1] : zero(T)),
         (s <= 4 ? (rho[j]*h)^(4-s)/fac[4-s+1] : zero(T)),
         # w coefficients: h^(nd-s+1) * rho[j]^(k+nd-s) / (k+nd-s)! for k = 1:nc
         h^(nd-s+1)*rho[j]^(1+nd-s)/fac[1+nd-s+1],
         h^(nd-s+1)*rho[j]^(2+nd-s)/fac[2+nd-s+1],
         h^(nd-s+1)*rho[j]^(3+nd-s)/fac[3+nd-s+1],
         h^(nd-s+1)*rho[j]^(4+nd-s)/fac[4+nd-s+1])
        for i in 1:nh, j in 1:nc, s in 1:nd
    ]

    # continuity coefficients: depend on s
    cont_itr = [
        (i, s,
         # v coefficients
         (s <= 1 ? h^(1-s)/fac[1-s+1] : zero(T)),
         (s <= 2 ? h^(2-s)/fac[2-s+1] : zero(T)),
         (s <= 3 ? h^(3-s)/fac[3-s+1] : zero(T)),
         (s <= 4 ? h^(4-s)/fac[4-s+1] : zero(T)),
         # w coefficients
         h^(nd-s+1)/fac[1+nd-s+1],
         h^(nd-s+1)/fac[2+nd-s+1],
         h^(nd-s+1)/fac[3+nd-s+1],
         h^(nd-s+1)/fac[4+nd-s+1])
        for i in 1:nh-1, s in 1:nd
    ]

    # collocation physics coefficients: depend on j
    coll_itr = [
        (i, j,
         rho[j]^0/fac[1], rho[j]^1/fac[2], rho[j]^2/fac[3], rho[j]^3/fac[4])
        for i in 1:nh, j in 1:nc
    ]

    # right BC coefficients
    bc3_cv = [h^(k-1)/fac[k] for k in 1:nd]
    bc3_cw = [h^nd/fac[k+nd] for k in 1:nc]
    bc4_cv = [h^(k-2)/fac[k-1] for k in 2:nd]
    bc4_cw = [h^(nd-1)/fac[k+nd-1] for k in 1:nc]

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@add_variable(core, v, nh, nd; start = v0)
    ExaModels.@add_variable(core, w, nh, nc; start = 0.0)
    ExaModels.@add_variable(core, uc, nh, nc, nd; start = uc0)
    ExaModels.@add_variable(core, Duc, nh, nc, nd; start = 0.0)

    # Constant objective
    ExaModels.@add_objective(core, one(T) for _i in 1:1)

    # con1: uc[i,j,s] - v[i,s] - h * sum(w[i,k] * rho[j]^k / k!)
    ExaModels.@add_constraint(core, c1,
        uc[i, j, s] - v[i, s] - (a1*w[i,1] + a2*w[i,2] + a3*w[i,3] + a4*w[i,4])
        for (i, j, s, a1, a2, a3, a4) in con1_itr
    )

    # con2: Duc[i,j,s] - sum(v[i,k]*Bv[k]) - sum(w[i,k]*Bw[k])
    ExaModels.@add_constraint(core, c2,
        Duc[i, j, s] - (bv1*v[i,1] + bv2*v[i,2] + bv3*v[i,3] + bv4*v[i,4])
                     - (bw1*w[i,1] + bw2*w[i,2] + bw3*w[i,3] + bw4*w[i,4])
        for (i, j, s, bv1, bv2, bv3, bv4, bw1, bw2, bw3, bw4) in con2_itr
    )

    # Boundary conditions
    ExaModels.@add_constraint(core, bc1, v[1, 1] - bc[1, 1])
    ExaModels.@add_constraint(core, bc2, v[1, 2] - bc[2, 1])

    ExaModels.@add_constraint(core, bc3,
        bc3_cv[1]*v[nh,1] + bc3_cv[2]*v[nh,2] + bc3_cv[3]*v[nh,3] + bc3_cv[4]*v[nh,4] +
        bc3_cw[1]*w[nh,1] + bc3_cw[2]*w[nh,2] + bc3_cw[3]*w[nh,3] + bc3_cw[4]*w[nh,4] - bc[1, 2]
    )

    ExaModels.@add_constraint(core, bc4,
        bc4_cv[1]*v[nh,2] + bc4_cv[2]*v[nh,3] + bc4_cv[3]*v[nh,4] +
        bc4_cw[1]*w[nh,1] + bc4_cw[2]*w[nh,2] + bc4_cw[3]*w[nh,3] + bc4_cw[4]*w[nh,4] - bc[2, 2]
    )

    # Continuity
    ExaModels.@add_constraint(core, c3,
        cv1*v[i,1] + cv2*v[i,2] + cv3*v[i,3] + cv4*v[i,4] +
        cw1*w[i,1] + cw2*w[i,2] + cw3*w[i,3] + cw4*w[i,4] - v[i+1, s]
        for (i, s, cv1, cv2, cv3, cv4, cw1, cw2, cw3, cw4) in cont_itr
    )

    # Collocation physics
    ExaModels.@add_constraint(core, c4,
        e1*w[i,1] + e2*w[i,2] + e3*w[i,3] + e4*w[i,4] -
        R * (Duc[i, j, 2] * Duc[i, j, 3] - Duc[i, j, 1] * Duc[i, j, 4])
        for (i, j, e1, e2, e3, e4) in coll_itr
    )

    return ExaModels.ExaModel(core; kwargs...)
end
