# Journal Bearing Problem

function bearing_model end

const bearing_b = 10    # grid is (0,2*pi)x(0,2*b)
const bearing_e = 0.1   # eccentricity

@inline bearing_wq(i, hx)   = (1.0 + bearing_e*cos((i-1)*hx))^3
@inline bearing_v0_ij(i, hx) = max(sin((i-1)*hx), 0.0)
