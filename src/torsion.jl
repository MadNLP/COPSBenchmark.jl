# Torsion Problem

function torsion_model end

const torsion_c_val = 5.0

function torsion_D(nx, ny; T = Float64)
    hx = T(1.0 / (nx + 1.0))
    hy = T(1.0 / (ny + 1.0))
    return T[min(min(i, nx-i+1)*hx, min(j, ny-j+1)*hy) for i in 0:nx+1, j in 0:ny+1]
end
