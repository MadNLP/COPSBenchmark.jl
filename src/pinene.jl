# Isomerization of Alpha-Pinene Problem
# Collocation formulation

function pinene_model end

const pinene_nc = 3
const pinene_ne = 5
const pinene_np = 5
const pinene_nm = 8
const pinene_rho = [0.11270166537926, 0.5, 0.88729833462074]
const pinene_bc  = [100.0, 0.0, 0.0, 0.0, 0.0]
const pinene_tau = [1230.0, 3060.0, 4920.0, 7800.0, 10680.0, 15030.0, 22620.0, 36420.0]
const pinene_z_data = [
    88.35    7.3     2.3     0.4     1.75;
    76.4    15.6     4.5     0.7     2.8;
    65.1    23.1     5.3     1.1     5.8;
    50.4    32.9     6.0     1.5     9.3;
    37.5    42.7     6.0     1.9    12.0;
    25.9    49.1     5.9     2.2    17.0;
    14.0    57.4     5.1     2.6    21.0;
    4.5    63.1     3.8     2.9    25.7;
]

function pinene_data(nh; T = Float64)
    tf   = pinene_tau[pinene_nm]
    h    = tf / nh
    t    = [(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, floor(pinene_tau[i]/h)+1) for i in 1:pinene_nm]
    z    = T.(pinene_z_data)
    v0   = zeros(T, nh, pinene_ne)
    for i in 1:itau[1], s in 1:pinene_ne
        v0[i, s] = pinene_bc[s]
    end
    for j in 2:pinene_nm, i in itau[j-1]+1:itau[j], s in 1:pinene_ne
        v0[i, s] = z[j, s]
    end
    for i in itau[pinene_nm]+1:nh, s in 1:pinene_ne
        v0[i, s] = z[pinene_nm, s]
    end
    return (; h, t, itau, z, v0)
end

@inline pinene_ode1(theta, uc1)            = -(theta[1]+theta[2])*uc1
@inline pinene_ode2(theta, uc1)            = theta[1]*uc1
@inline pinene_ode3(theta, uc1, uc3, uc5)  = theta[2]*uc1 - (theta[3]+theta[4])*uc3 + theta[5]*uc5
@inline pinene_ode4(theta, uc3)            = theta[3]*uc3
@inline pinene_ode5(theta, uc3, uc5)       = theta[4]*uc3 - theta[5]*uc5
