# Electrons on a Sphere Problem

function elec_model end

@inline elec_coulomb_potential(xi, yi, zi, xj, yj, zj) =
    1.0 / sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
