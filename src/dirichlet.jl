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

function dirichlet_model(nh)
    dom = PDEDiscretizationDomain(nh, CIRCLE_DOMAIN)
    pb = PDEProblem(
        0.01,
        fill(1.0, dom.NODES),
        fill(1.0, dom.NODES),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return _transition_state_model(pb, dom)
end

