# Catalyst Mixing Problem
# Collocation formulation

function catmix_model end

const catmix_ne = 2
const catmix_nc = 3
const catmix_tf = 1
const catmix_rho = [0.11270166537926, 0.50000000000000, 0.88729833462074]
const catmix_bc  = [1.0, 0.0]
const catmix_alpha = 0.0  # smoothing parameter

@inline catmix_ode1(u, pp1, pp2) = u * (10.0*pp2 - pp1)
@inline catmix_ode2(u, pp1, pp2) = u * (pp1 - 10.0*pp2) - (1 - u)*pp2
