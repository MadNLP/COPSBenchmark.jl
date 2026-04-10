# Minimal Surface with Obstacle Problem

function minsurf_model end

@inline minsurf_v0(x_mesh, i, j) = 1 - (2 * x_mesh[i] - 1)^2
