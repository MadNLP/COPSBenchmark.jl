# Catalytic Cracking of Gas Oil Problem
# Collocation formulation

function gasoil_model end

const gasoil_nc = 4
const gasoil_ne = 2
const gasoil_np = 3
const gasoil_nm = 21
const gasoil_rho = [0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
const gasoil_bc  = [1, 1, 2, 0]  # ODE initial conditions (only first ne used for v0)
const gasoil_tau = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.150, 0.175, 0.20, 0.225, 0.250,
                    0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65, 0.75, 0.85, 0.95]
const gasoil_z_data = [
    1.0000        0;
    0.8105   0.2000;
    0.6208   0.2886;
    0.5258   0.3010;
    0.4345   0.3215;
    0.3903   0.3123;
    0.3342   0.2716;
    0.3034   0.2551;
    0.2735   0.2258;
    0.2405   0.1959;
    0.2283   0.1789;
    0.2071   0.1457;
    0.1669   0.1198;
    0.1530   0.0909;
    0.1339   0.0719;
    0.1265   0.0561;
    0.1200   0.0460;
    0.0990   0.0280;
    0.0870   0.0190;
    0.0770   0.0140;
    0.0690   0.0100;
]

function gasoil_data(nh; T = Float64)
    tf   = gasoil_tau[gasoil_nm]
    h    = tf / nh
    t    = [(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, floor(gasoil_tau[i]/h)+1) for i in 1:gasoil_nm]
    z    = T.(gasoil_z_data)
    v0   = zeros(T, nh, gasoil_ne)
    for i in 1:itau[1], s in 1:gasoil_ne
        v0[i, s] = gasoil_bc[s]
    end
    for j in 2:gasoil_nm, i in itau[j-1]+1:itau[j], s in 1:gasoil_ne
        v0[i, s] = z[j, s]
    end
    for i in itau[gasoil_nm]+1:nh, s in 1:gasoil_ne
        v0[i, s] = z[gasoil_nm, s]
    end
    return (; h, t, itau, z, v0)
end

@inline gasoil_ode1(theta, uc1)       = -(theta[1]+theta[3])*uc1^2
@inline gasoil_ode2(theta, uc1, uc2)  = theta[1]*uc1^2 - theta[2]*uc2
