# Hanging Chain Problem

function chain_model end

const chain_L    = 4
const chain_a    = 1
const chain_b    = 3
const chain_tmin = chain_b > chain_a ? 1/4 : 3/4
const chain_tf   = 1.0

@inline chain_nh(n)             = max(2, div(n - 4, 4))
@inline chain_u_start(k, nh)   = 4 * abs(chain_b - chain_a) * (k / nh - chain_tmin)
@inline chain_x1_start(k, nh)  = 4 * abs(chain_b - chain_a) * k / nh * (1/2 * k / nh - chain_tmin) + chain_a
@inline chain_x2_start(k, nh)  = chain_x1_start(k, nh) * chain_u_start(k, nh)
