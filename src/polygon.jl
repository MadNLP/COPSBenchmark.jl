# Polygon of Maximal Area Problem

function polygon_model end

@inline polygon_N(n)                   = div(n, 2)
@inline polygon_area_term(ri, ri1, ti, ti1) = -0.5 * ri * ri1 * sin(ti1 - ti)
