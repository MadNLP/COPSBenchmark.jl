# Flow in a Channel Problem
# Collocation formulation

function channel_model end

const channel_nc = 4
const channel_nd = 4
const channel_R  = 10.0
const channel_tf = 1.0
const channel_bc = [0.0 1.0; 0.0 0.0]
const channel_rho = [0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]

function channel_v0(nh; T = Float64)
    h  = T(channel_tf) / nh
    t  = T[(i-1)*h for i in 1:nh+1]
    v0 = zeros(T, nh, channel_nd)
    for i in 1:nh
        v0[i, 1] = t[i]^2*(3 - 2*t[i])
        v0[i, 2] = 6*t[i]*(1 - t[i])
        v0[i, 3] = 6*(1 - 2*t[i])
        v0[i, 4] = -12
    end
    return (; h, t, v0)
end

@inline channel_collocation_rhs(R, Duc1, Duc2, Duc3, Duc4) = R * (Duc2 * Duc3 - Duc1 * Duc4)
