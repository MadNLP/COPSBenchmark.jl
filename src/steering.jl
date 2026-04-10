# Rocket Steering Problem

function steering_model end

const steering_a     = 100.0
const steering_u_min = -pi/2.0
const steering_u_max =  pi/2.0
const steering_xs    = zeros(4)
const steering_xf    = [NaN, 5.0, 45.0, 0.0]

@inline function steering_x0(k, i, nh)
    if i == 1 || i == 4
        return 0.0
    elseif i == 2
        return 5*k/nh
    elseif i == 3
        return 45.0*k/nh
    else
        return 0.0
    end
end
