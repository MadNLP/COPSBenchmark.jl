# Methanol-to-Hydrocarbons Problem
# Collocation formulation

function methanol_model end

const methanol_ne = 3
const methanol_np = 5
const methanol_nc = 3
const methanol_nm = 17
const methanol_rho = [0.11270166537926, 0.5, 0.88729833462074]
const methanol_bc  = [1.0, 0.0, 0.0]
const methanol_tau = [
    0.,
    0.050,
    0.065,
    0.080,
    0.123,
    0.233,
    0.273,
    0.354,
    0.397,
    0.418,
    0.502,
    0.553,
    0.681,
    0.750,
    0.916,
    0.937,
    1.122,
]
const methanol_z_data = [
    1.0000  0.0000  0.0000;
    0.7085  0.1621  0.0811;
    0.5971  0.1855  0.0965;
    0.5537  0.1989  0.1198;
    0.3684  0.2845  0.1535;
    0.1712  0.3491  0.2097;
    0.1198  0.3098  0.2628;
    0.0747  0.3576  0.2467;
    0.0529  0.3347  0.2884;
    0.0415  0.3388  0.2757;
    0.0261  0.3557  0.3167;
    0.0208  0.3483  0.2954;
    0.0085  0.3836  0.2950;
    0.0053  0.3611  0.2937;
    0.0019  0.3609  0.2831;
    0.0018  0.3485  0.2846;
    0.0006  0.3698  0.2899;
]

function methanol_data(nh; T = Float64)
    tf   = methanol_tau[methanol_nm]
    h    = tf / nh
    t    = [(i-1)*h for i in 1:nh+1]
    itau = Int[min(nh, floor(methanol_tau[i]/h)+1) for i in 1:methanol_nm]
    z    = T.(methanol_z_data)
    v0   = zeros(T, nh, methanol_ne)
    for i in 1:itau[1], s in 1:methanol_ne
        v0[i, s] = methanol_bc[s]
    end
    for j in 2:methanol_nm, i in itau[j-1]+1:itau[j], s in 1:methanol_ne
        v0[i, s] = z[j, s]
    end
    for i in itau[methanol_nm]+1:nh, s in 1:methanol_ne
        v0[i, s] = z[methanol_nm, s]
    end
    v0 .= 0.001
    return (; h, t, itau, z, v0)
end

@inline methanol_ode1(theta, uc1, uc2) = -(2*theta[2] - (theta[1]*uc2)/((theta[2]+theta[5])*uc1+uc2) + theta[3] + theta[4])*uc1
@inline methanol_ode2(theta, uc1, uc2) = (theta[1]*uc1*(theta[2]*uc1-uc2))/((theta[2]+theta[5])*uc1+uc2) + theta[3]*uc1
@inline methanol_ode3(theta, uc1, uc2) = (theta[1]*uc1*(uc2+theta[5]*uc1))/((theta[2]+theta[5])*uc1+uc2) + theta[4]*uc1
