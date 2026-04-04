# Flow in a Channel Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function COPSBenchmark.channel_model(::ExaModelsBackend, nh; T = Float64, backend = nothing, kwargs...)
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
    fac = [factorial(k) for k in 0:nc+nd]

    core = ExaModels.ExaCore(T; backend = backend)

    ExaModels.@var(core, v, nh, nd; start = v0)
    ExaModels.@var(core, w, nh, nc; start = 0.0)
    ExaModels.@var(core, uc, nh, nc, nd; start = uc0)
    ExaModels.@var(core, Duc, nh, nc, nd; start = 0.0)

    # con1: uc[i,j,s] = v[i,s] + h * sum_{k=1}^{nc}(w[i,k] * rho[j]^k / k!)
    # Expand over j so rho[j] is a concrete constant per pattern
    for j_val in 1:nc
        rho_j = rho[j_val]
        A = [rho_j^k / fac[k+1] for k in 1:nc]  # A[k] = rho_j^k / k!
        ExaModels.@con(
            core,
            uc[i, j_val, s] - v[i, s] - h * (w[i,1]*A[1] + w[i,2]*A[2] + w[i,3]*A[3] + w[i,4]*A[4])
            for i in 1:nh, s in 1:nd
        )
    end

    # con2: Duc[i,j,s] = sum_{k=s}^{nd}(v[i,k]*(rho[j]*h)^(k-s)/(k-s)!)
    #                  + h^(nd-s+1) * sum_{k=1}^{nc}(w[i,k]*rho[j]^(k+nd-s)/(k+nd-s)!)
    # Expand over j and s so all coefficients are concrete
    for j_val in 1:nc, s_val in 1:nd
        rho_j = rho[j_val]
        # v_coefs[m] for k = s_val+m-1, m in 1:nd-s_val+1
        Bv = [(rho_j * h)^(k - s_val) / fac[k - s_val + 1] for k in s_val:nd]
        # w_coefs[k] for k in 1:nc
        Bw = [h^(nd - s_val + 1) * rho_j^(k + nd - s_val) / fac[k + nd - s_val + 1] for k in 1:nc]
        ExaModels.@con(
            core,
            Duc[i, j_val, s_val]
            - sum(v[i, k] * Bv[k - s_val + 1] for k in s_val:nd)
            - (Bw[1]*w[i,1] + Bw[2]*w[i,2] + Bw[3]*w[i,3] + Bw[4]*w[i,4])
            for i in 1:nh
        )
    end

    # Boundary conditions
    ExaModels.@con(core, bc1, v[1, 1] - bc[1, 1])
    ExaModels.@con(core, bc2, v[1, 2] - bc[2, 1])

    # Right boundary: sum_{k=1}^{nd}(v[nh,k]*h^(k-1)/(k-1)!) + h^nd*sum_{k=1}^{nc}(w[nh,k]/(k+nd-1)!) = bc[1,2]
    let Cv = [h^(k-1) / fac[k] for k in 1:nd],
        Cw = [h^nd / fac[k+nd] for k in 1:nc]
        ExaModels.@con(
            core, bc3,
            sum(v[nh, k] * Cv[k] for k in 1:nd) +
            (Cw[1]*w[nh,1] + Cw[2]*w[nh,2] + Cw[3]*w[nh,3] + Cw[4]*w[nh,4]) - bc[1, 2]
        )
    end

    # Right boundary: sum_{k=2}^{nd}(v[nh,k]*h^(k-2)/(k-2)!) + h^(nd-1)*sum_{k=1}^{nc}(w[nh,k]/(k+nd-2)!) = bc[2,2]
    let Cv = [h^(k-2) / fac[k-1] for k in 2:nd],
        Cw = [h^(nd-1) / fac[k+nd-1] for k in 1:nc]
        ExaModels.@con(
            core, bc4,
            sum(v[nh, k] * Cv[k-1] for k in 2:nd) +
            (Cw[1]*w[nh,1] + Cw[2]*w[nh,2] + Cw[3]*w[nh,3] + Cw[4]*w[nh,4]) - bc[2, 2]
        )
    end

    # Continuity: sum_{k=s}^{nd}(v[i,k]*h^(k-s)/(k-s)!) + h^(nd-s+1)*sum_{k=1}^{nc}(w[i,k]/(k+nd-s)!) = v[i+1,s]
    # Expand over s so coefficients are concrete
    for s_val in 1:nd
        Dv = [h^(k - s_val) / fac[k - s_val + 1] for k in s_val:nd]
        Dw = [h^(nd - s_val + 1) / fac[k + nd - s_val + 1] for k in 1:nc]
        ExaModels.@con(
            core,
            sum(v[i, k] * Dv[k - s_val + 1] for k in s_val:nd) +
            (Dw[1]*w[i,1] + Dw[2]*w[i,2] + Dw[3]*w[i,3] + Dw[4]*w[i,4]) - v[i+1, s_val]
            for i in 1:nh-1
        )
    end

    # Collocation physics: sum_{k=1}^{nc}(w[i,k]*rho[j]^(k-1)/(k-1)!) = R*(Duc[i,j,2]*Duc[i,j,3] - Duc[i,j,1]*Duc[i,j,4])
    # Expand over j so rho[j] is concrete
    for j_val in 1:nc
        rho_j = rho[j_val]
        E = [rho_j^(k-1) / fac[k] for k in 1:nc]  # E[k] = rho_j^(k-1) / (k-1)!
        ExaModels.@con(
            core,
            (E[1]*w[i,1] + E[2]*w[i,2] + E[3]*w[i,3] + E[4]*w[i,4]) -
            R * (Duc[i, j_val, 2] * Duc[i, j_val, 3] - Duc[i, j_val, 1] * Duc[i, j_val, 4])
            for i in 1:nh
        )
    end

    return ExaModels.ExaModel(core; kwargs...)
end
