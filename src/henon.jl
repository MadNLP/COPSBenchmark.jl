# Todd Munson
# COPS 3.1 - March 2004

###############################################################################
# Generalized Lane-Emden-Fowler-Henon equation with homogeneous Dirichlet
# boundary conditions:
#
#   a*Laplacian - b(x)*u + c(x)*u - d(x) = 0	for x in  int Omega
#   u(x) = 0					for x in bdry Omega
#
# Equation classes considered:
#   Lane-Emden-Fowler: c(x) = 1
#   Henon            : c(x) = |x|^(2*l)
###############################################################################

function henon_model(nh)
    dom = PDEDiscretizationDomain(nh, CIRCLE_REC_DOMAIN)
    pb = PDEProblem(
        1.0,
        fill(0.0, dom.NODES),
        sqrt.(dom.COORDS[:, 1].^2 .+ dom.COORDS[:, 2].^2),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return _transition_state_model(pb, dom)
end

