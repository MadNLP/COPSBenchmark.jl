# Cam Shape Problem

function camshape_model end

const camshape_R_v   = 1.0
const camshape_R_max = 2.0
const camshape_R_min = 1.0
const camshape_alpha = 1.5

@inline camshape_d_theta(n)           = 2*pi/(5*(n+1))
@inline camshape_convexity(r1, r2, r3, d_theta) = -r1*r2 - r2*r3 + 2*r1*r3*cos(d_theta)
