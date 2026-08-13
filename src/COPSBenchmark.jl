module COPSBenchmark

using Random

abstract type AbstractModelerBackend end
struct JuMPBackend <: AbstractModelerBackend end
struct ExaModelsBackend <: AbstractModelerBackend end

# COPS Instances
function bearing_model end
function bearing_recipe end
function bearing_args end
# Two sizes, so a two-argument deferred call; the recipe then carries only the
# two integers, which is what a builder ABI can accept.
function bearing_data(nx, ny)
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
    return (; v0, lower, upper, lin)
end
function camshape_model end
function camshape_recipe end
function camshape_args end
function catmix_model end
function catmix_recipe end
function catmix_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function catmix_data(nh; T = Float64)
    ne = 2
    nc = 3
    # Concrete element type, not the keyword `T`: inside a compiled library `T`
    # is a runtime value and `T[...]` leaves `_array_for(T::Type, ...)`
    # unresolved for `--trim=safe`.  The core converts on append, so the values
    # land in the model's own type either way.
    v_start  = Float64[mod(j, ne) for i in 1:nh, j in 1:ne]
    pp_start = Float64[mod(k, ne) for i in 1:nh, j in 1:nc, k in 1:ne]
    return (; v_start, pp_start)
end
function chain_model end
function chain_recipe end
function chain_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function chain_data(nh; T = Float64)
    a = 1
    b = 3
    tmin = b > a ? 1 / 4 : 3 / 4
    u0  = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1]
    x10 = [4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a for k in 1:nh+1]
    x20 = [(4 * abs(b - a) * k / nh * (1 / 2 * k / nh - tmin) + a) *
           (4 * abs(b - a) * (k / nh - tmin)) for k in 1:nh+1]
    x30 = [4 * abs(b - a) * (k / nh - tmin) for k in 1:nh+1]
    return (; u0, x10, x20, x30)
end
function channel_model end
function channel_recipe end
function channel_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function channel_data(nh; T = Float64)
    nc = 4
    nd = 4
    tf = T(1.0)
    h = tf / nh

    rho = Float64[0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    t = Float64[(i-1)*h for i in 1:nh+1]

    # Initial value
    v0 = zeros(T, nh, nd)
    for i in 1:nh
        v0[i, 1] = t[i]^2*(3 - 2*t[i])
        v0[i, 2] = 6*t[i]*(1 - t[i])
        v0[i, 3] = 6*(1 - 2*t[i])
        v0[i, 4] = -12
    end
    uc0 = Float64[v0[i, s] for i in 1:nh, j in 1:nc, s in 1:nd]

    # fac[k+1] = k!
    fac = Float64[factorial(k) for k in 0:nc+nd]

    con1_itr = Tuple{Int,Int,Int,T,T,T,T}[
        (i, j, s, h*rho[j]^1/fac[2], h*rho[j]^2/fac[3], h*rho[j]^3/fac[4], h*rho[j]^4/fac[5])
        for i in 1:nh, j in 1:nc, s in 1:nd
    ]

    con2_itr = Tuple{Int,Int,Int,T,T,T,T,T,T,T,T}[
        (i, j, s,
         (s <= 1 ? (rho[j]*h)^(1-s)/fac[1-s+1] : zero(T)),
         (s <= 2 ? (rho[j]*h)^(2-s)/fac[2-s+1] : zero(T)),
         (s <= 3 ? (rho[j]*h)^(3-s)/fac[3-s+1] : zero(T)),
         (s <= 4 ? (rho[j]*h)^(4-s)/fac[4-s+1] : zero(T)),
         h^(nd-s+1)*rho[j]^(1+nd-s)/fac[1+nd-s+1],
         h^(nd-s+1)*rho[j]^(2+nd-s)/fac[2+nd-s+1],
         h^(nd-s+1)*rho[j]^(3+nd-s)/fac[3+nd-s+1],
         h^(nd-s+1)*rho[j]^(4+nd-s)/fac[4+nd-s+1])
        for i in 1:nh, j in 1:nc, s in 1:nd
    ]

    cont_itr = Tuple{Int,Int,T,T,T,T,T,T,T,T}[
        (i, s,
         (s <= 1 ? h^(1-s)/fac[1-s+1] : zero(T)),
         (s <= 2 ? h^(2-s)/fac[2-s+1] : zero(T)),
         (s <= 3 ? h^(3-s)/fac[3-s+1] : zero(T)),
         (s <= 4 ? h^(4-s)/fac[4-s+1] : zero(T)),
         h^(nd-s+1)/fac[1+nd-s+1],
         h^(nd-s+1)/fac[2+nd-s+1],
         h^(nd-s+1)/fac[3+nd-s+1],
         h^(nd-s+1)/fac[4+nd-s+1])
        for i in 1:nh-1, s in 1:nd
    ]

    coll_itr = Tuple{Int,Int,T,T,T,T}[
        (i, j, rho[j]^0/fac[1], rho[j]^1/fac[2], rho[j]^2/fac[3], rho[j]^3/fac[4])
        for i in 1:nh, j in 1:nc
    ]

    # right BC coefficients, carried with the row index they apply to
    bc3_cv = [h^(k-1)/fac[k] for k in 1:nd]
    bc3_cw = [h^nd/fac[k+nd] for k in 1:nc]
    bc4_cv = [h^(k-2)/fac[k-1] for k in 2:nd]
    bc4_cw = [h^(nd-1)/fac[k+nd-1] for k in 1:nc]
    bc3_itr = Tuple{Int,T,T,T,T,T,T,T,T}[(nh, bc3_cv[1], bc3_cv[2], bc3_cv[3], bc3_cv[4],
                    bc3_cw[1], bc3_cw[2], bc3_cw[3], bc3_cw[4])]
    bc4_itr = Tuple{Int,T,T,T,T,T,T,T}[(nh, bc4_cv[1], bc4_cv[2], bc4_cv[3],
                    bc4_cw[1], bc4_cw[2], bc4_cw[3], bc4_cw[4])]
    return (; v0, uc0, con1_itr, con2_itr, cont_itr, coll_itr, bc3_itr, bc4_itr)
end
function clnlbeam_model end
function dirichlet_model end
function elec_model end
function elec_recipe end
function elec_args end
# The draw is data, and the seed is part of it.  `Fix2` on a NAMED function keeps
# every type in the core named, and bakes the seed into a compiled library --
# which is what a library that instantiates from one integer must do.
function elec_data(np, seed)
    Random.seed!(seed)
    theta = (2pi) .* rand(np)
    phi = pi .* rand(np)
    x0 = Float64[cos(theta[i])*sin(phi[i]) for i=1:np]
    y0 = Float64[sin(theta[i])*sin(phi[i]) for i=1:np]
    z0 = Float64[cos(phi[i]) for i=1:np]
    itr = Tuple{Int,Int}[(i,j) for i in 1:np-1 for j in i+1:np]
    return (; x0, y0, z0, itr)
end
function gasoil_model end
function gasoil_recipe end
function gasoil_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function gasoil_data(nh; T = Float64)
    nc = 4
    ne = 2
    nm = 21
    rho = Float64[0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    bc = [1, 1, 2, 0]
    tau = Float64[0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.150, 0.175, 0.20, 0.225, 0.250, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85, 0.95]
    tf = tau[nm]
    h = tf / T(nh)
    t = Float64[T(i-1)*h for i in 1:nh+1]

    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(Float64[
        1.0000, 0.0000,
        0.8105, 0.2000,
        0.6208, 0.2886,
        0.5258, 0.3010,
        0.4345, 0.3215,
        0.3903, 0.3123,
        0.3342, 0.2716,
        0.3034, 0.2551,
        0.2735, 0.2258,
        0.2405, 0.1959,
        0.2283, 0.1789,
        0.2071, 0.1457,
        0.1669, 0.1198,
        0.1530, 0.0909,
        0.1339, 0.0719,
        0.1265, 0.0561,
        0.1200, 0.0460,
        0.0990, 0.0280,
        0.0870, 0.0190,
        0.0770, 0.0140,
        0.0690, 0.0100,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = T(bc[s])
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start  = [v0[i, s] for i = 1:nh, s = 1:ne]
    uc_start = [v0[i, s] for i = 1:nh, j = 1:nc, s = 1:ne]
    itr = Tuple{Int,Int,Int,T,T,T}[(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    return (; v_start, uc_start, itr)
end
function glider_model end
function glider_recipe end
function glider_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function glider_data(nh; T = Float64)
    x_0 = T(0); y_0 = T(1000); y_f = T(900); vx_0 = T(13.23)
    inv_nh = T(1) / T(nh)
    x_start = [x_0 + vx_0*(T(k)*inv_nh) for k in 0:nh]
    y_start = [y_0 + (T(k)*inv_nh)*(y_f - y_0) for k in 0:nh]
    return (; x_start, y_start)
end
function henon_model end
function lane_emden_model end
function marine_model end
function marine_recipe end
function marine_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function marine_data(nh; T = Float64)
    ne = 8
    nm = 21
    tau = collect(T, range(T(0), T(10), 21))
    tf = tau[nm]
    h = tf / T(nh)
    t = Float64[T(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(Float64[
        20000.0, 17000.0, 10000.0, 15000.0, 12000.0,  9000.0,  7000.0,  3000.0,
        12445.0, 15411.0, 13040.0, 13338.0, 13484.0,  8426.0,  6615.0,  4022.0,
        7705.0, 13074.0, 14623.0, 11976.0, 12453.0,  9272.0,  6891.0,  5020.0,
        4664.0,  8579.0, 12434.0, 12603.0, 11738.0,  9710.0,  6821.0,  5722.0,
        2977.0,  7053.0, 11219.0, 11340.0, 13665.0,  8534.0,  6242.0,  5695.0,
        1769.0,  5054.0, 10065.0, 11232.0, 12112.0,  9600.0,  6647.0,  7034.0,
        943.0,  3907.0,  9473.0, 10334.0, 11115.0,  8826.0,  6842.0,  7348.0,
        581.0,  2624.0,  7421.0, 10297.0, 12427.0,  8747.0,  7199.0,  7684.0,
        355.0,  1744.0,  5369.0,  7748.0, 10057.0,  8698.0,  6542.0,  7410.0,
        223.0,  1272.0,  4713.0,  6869.0,  9564.0,  8766.0,  6810.0,  6961.0,
        137.0,   821.0,  3451.0,  6050.0,  8671.0,  8291.0,  6827.0,  7525.0,
        87.0,   577.0,  2649.0,  5454.0,  8430.0,  7411.0,  6423.0,  8388.0,
        49.0,   337.0,  2058.0,  4115.0,  7435.0,  7627.0,  6268.0,  7189.0,
        32.0,   228.0,  1440.0,  3790.0,  6474.0,  6658.0,  5859.0,  7467.0,
        17.0,   168.0,  1178.0,  3087.0,  6524.0,  5880.0,  5562.0,  7144.0,
        11.0,    99.0,   919.0,  2596.0,  5360.0,  5762.0,  4480.0,  7256.0,
        7.0,    65.0,   647.0,  1873.0,  4556.0,  5058.0,  4944.0,  7538.0,
        4.0,    44.0,   509.0,  1571.0,  4009.0,  4527.0,  4233.0,  6649.0,
        2.0,    27.0,   345.0,  1227.0,  3677.0,  4229.0,  3805.0,  6378.0,
        1.0,    20.0,   231.0,   934.0,  3197.0,  3695.0,  3159.0,  6454.0,
        1.0,    12.0,   198.0,   707.0,  2562.0,  3163.0,  3232.0,  5566.0,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = z[1, s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start = [v0[i, s] for i=1:nh, s=1:ne]
    itr = Tuple{Int,Int,Int,T,T,T}[(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    return (; v_start, itr)
end
function methanol_model end
function methanol_recipe end
function methanol_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function methanol_data(nh; T = Float64)
    ne = 3
    nm = 17
    tau = Float64[
        0., 0.050, 0.065, 0.080, 0.123, 0.233, 0.273, 0.354, 0.397, 0.418,
        0.502, 0.553, 0.681, 0.750, 0.916, 0.937, 1.122,
    ]
    tf = tau[nm]
    h = tf / T(nh)
    t = Float64[T(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(Float64[
        1.0000, 0.0000, 0.0000,
        0.7085, 0.1621, 0.0811,
        0.5971, 0.1855, 0.0965,
        0.5537, 0.1989, 0.1198,
        0.3684, 0.2845, 0.1535,
        0.1712, 0.3491, 0.2097,
        0.1198, 0.3098, 0.2628,
        0.0747, 0.3576, 0.2467,
        0.0529, 0.3347, 0.2884,
        0.0415, 0.3388, 0.2757,
        0.0261, 0.3557, 0.3167,
        0.0208, 0.3483, 0.2954,
        0.0085, 0.3836, 0.2950,
        0.0053, 0.3611, 0.2937,
        0.0019, 0.3609, 0.2831,
        0.0018, 0.3485, 0.2846,
        0.0006, 0.3698, 0.2899,
    ], ne, nm)'

    con1_matrix = Tuple{Int,Int,Int,T,T,T}[(j, s, itau[j], tau[j], z[j,s], t[itau[j]]) for j in 1:nm, s in 1:ne]
    return (; con1_matrix)
end
function minsurf_model end
function minsurf_recipe end
function minsurf_args end
# Two sizes, so a two-argument deferred call; the recipe then carries only the
# two integers, which is what a builder ABI can accept.
function minsurf_data(nx, ny)
    x_mesh = LinRange(0, 1, nx + 2) # coordinates of the mesh points x

    v0 = zeros(nx + 2, ny + 2) # Surface matrix initialization
    for i = 1:(nx + 2), j = 1:(ny + 2)
        v0[i, j] = 1 - (2 * x_mesh[i] - 1)^2
    end

    hx = 1 / (nx + 1)
    hy = 1 / (ny + 1)
    c6_i = Int(floor(0.25 / hx)):Int(ceil(0.75 / hx))
    c6_j = Int(floor(0.25 / hy)):Int(ceil(0.75 / hy))
    return (; v0, c6_i, c6_j)
end
function pinene_model end
function pinene_recipe end
function pinene_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function pinene_data(nh; T = Float64)
    nc = 3
    ne = 5
    nm = 8
    bc = Float64[100, 0, 0, 0, 0]
    tau = Float64[1230, 3060, 4920, 7800, 10680, 15030, 22620, 36420]
    tf = tau[nm]
    h = tf / T(nh)
    t = Float64[T(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, Int(floor(tau[i]/h))+1) for i in 1:nm]

    z = reshape(Float64[
        88.35,  7.3,  2.3,  0.4,  1.75,
        76.4,  15.6,  4.5,  0.7,  2.8,
        65.1,  23.1,  5.3,  1.1,  5.8,
        50.4,  32.9,  6.0,  1.5,  9.3,
        37.5,  42.7,  6.0,  1.9, 12.0,
        25.9,  49.1,  5.9,  2.2, 17.0,
        14.0,  57.4,  5.1,  2.6, 21.0,
        4.5,  63.1,  3.8,  2.9, 25.7,
    ], ne, nm)'

    v0 = zeros(T, nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = bc[s]
    end
    for j in 2:nm, i =itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end

    v_start  = [v0[i, s] for i=1:nh, s=1:ne]
    uc_start = [v0[i, s] for i=1:nh, j=1:nc, s=1:ne]
    itr = Tuple{Int,Int,Int,T,T,T}[(j, s, itau[j], tau[j], t[itau[j]], z[j,s]) for j=1:nm, s in 1:ne]
    return (; v_start, uc_start, itr)
end
function polygon_model end
function polygon_recipe end
function polygon_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function polygon_data(N; T = Float64)
    θ0 = [i * π / (N - 1) - π / (N - 1) for i in 1:N]
    pairs = Tuple{Int,Int}[(i, j) for i in 1:N-1 for j in i+1:N]
    return (; θ0, pairs)
end
function robot_model end
function robot_recipe end
function robot_args end
# Derived tables, reached from the recipe through a deferred call so the
# recipe keeps a single placeholder -- which is all `P_new` carries.
function robot_data(nh; T = Float64)
    pi_T = T(pi)
    inv_nh = T(1) / T(nh)
    two_pi_3 = T(2) * pi_T / T(3)
    four_pi_3 = T(4) * pi_T / T(3)
    the_start     = [two_pi_3*(T(k)*inv_nh)^2 for k=1:nh+1]
    the_dot_start = [four_pi_3*(T(k)*inv_nh) for k=1:nh+1]
    return (; the_start, the_dot_start)
end
function rocket_model end
function rocket_recipe end
function rocket_args end
# Derived tables live here, not in `_args`: a recipe that reaches them through
# a deferred call keeps a single placeholder, which is all `P_new` carries.
function rocket_data(nh; T = Float64)
    m_0 = T(1)
    m_f = T(0.6) * m_0
    inv_nh = T(1) / T(nh)
    v_start = Float64[T(i)*inv_nh*(T(1) - T(i)*inv_nh) for i=0:nh]
    m_start = Float64[(m_f - m_0)*(T(i)*inv_nh) + m_0 for i=0:nh]
    return (; v_start, m_start)
end
function steering_model end
function steering_recipe end
function steering_args end
# Derived tables, reached through a deferred call so the recipe keeps a single
# placeholder.  Concrete element type: inside a compiled library the keyword `T`
# is a runtime value and `T[...]` will not trim.
function steering_data(nh)
    inv_nh = 1.0 / nh
    gen_x0(k, i) = i == 2 ? 5.0*k*inv_nh : (i == 3 ? 45.0*k*inv_nh : 0.0)
    x_start = Float64[gen_x0(i, j) for i = 1:nh+1, j = 1:4]
    return (; x_start)
end
function tetra_duct12_model end
function tetra_duct15_model end
function tetra_duct20_model end
function tetra_foam5_model end
function tetra_gear_model end
function tetra_hook_model end
function torsion_model end
function torsion_recipe end
function torsion_args end
# Two sizes, so a two-argument deferred call; the recipe then carries only the
# two integers, which is what a builder ABI can accept.
function torsion_data(nx, ny)
    hx = Float64(1.0 / (nx + 1.0))
    hy = Float64(1.0 / (ny + 1.0))
    D = Float64[min(min(i, nx-i+1)*hx, min(j, ny-j+1)*hy) for i in 0:nx+1, j in 0:ny+1]
    D_flat = [(k1, k2, D[k1, k2]) for k1 in 1:nx+2, k2 in 1:ny+2]
    lcon = [-d for (_, _, d) in D_flat]
    ucon = [d for (_, _, d) in D_flat]
    return (; D, D_flat, lcon, ucon)
end
function triangle_deer_model end
function triangle_pacman_model end
function triangle_turtle_model end

function transition_state_model end

include("pde_utils.jl")

end # module COPSBenchmark
