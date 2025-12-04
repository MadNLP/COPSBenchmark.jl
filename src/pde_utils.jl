
include("shapes/circle.jl")
include("shapes/circle_rec.jl")
include("shapes/rectangle.jl")

struct PDEProblem
    a::Float64
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    p::Vector{Float64}
end

struct PDEDiscretizationDomain
    NODES::Int
    ELEM::Int
    DIMEN::Int
    BREAK::Int
    AREA::Vector{Float64}
    TRIANG::Matrix{Int}
    COORDS::Matrix{Float64}
    BNDRY::Vector{Int}
    EDGE::Array{Float64, 3}
    US::Vector{Float64}
    UE::Vector{Float64}
end

function PDEDiscretizationDomain(nh, domain::Dict)
    NODES = domain[:NODES]
    ELEMS = domain[:ELEMS]
    COORDS = domain[:COORDS]
    ELEMS = domain[:ELEMS]
    BNDRY = domain[:BNDRY]
    DIMEN = 2
    BREAK = nh

    # Description of triangular elements
    TRIANG = domain[:TRIANG]
    # Edge lengths
    EDGE = [
        COORDS[TRIANG[e, mod(d1, DIMEN+1)+1], d2] - COORDS[TRIANG[e, d1], d2]
        for e in 1:ELEMS, d1 in 1:DIMEN+1, d2 in 1:DIMEN
    ]
    # Area of element
    AREA = [(EDGE[e, 1, 1]*EDGE[e, 2, 2] - EDGE[e, 1, 2]*EDGE[e, 2, 1]) / 2.0 for e in 1:ELEMS]
    US = domain[:US] # starting point
    UE = domain[:UE] # ending point

    return PDEDiscretizationDomain(
        NODES, ELEMS, DIMEN, BREAK,
        AREA, TRIANG, COORDS, BNDRY, EDGE, US, UE,
    )
end

function _update_values!(integral, energy, u, problem::PDEProblem, dom::PDEDiscretizationDomain)
    # Unpack values
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    TRIANG, EDGE = dom.TRIANG, dom.EDGE

    for b1 in 1:dom.BREAK+2, e1 in 1:dom.ELEM
        integral[b1, e1] =
            dom.AREA[e1]*(
            1 / (dom.DIMEN+1) *
                (sum((b[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]^2/2-
                            c[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]^(p[TRIANG[e1,c1]]+1)/(p[TRIANG[e1,c1]]+1)+
                            d[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]) for c1 in 1:dom.DIMEN+1)) +
            a / (8*dom.AREA[e1]^2)*(
                u[b1,TRIANG[e1,1]]^2*(EDGE[e1,2,1]^2 + EDGE[e1,2,2]^2) +
                u[b1,TRIANG[e1,2]]^2*(EDGE[e1,3,1]^2 + EDGE[e1,3,2]^2) +
                u[b1,TRIANG[e1,3]]^2*(EDGE[e1,1,1]^2 + EDGE[e1,1,2]^2) +
                2*u[b1,TRIANG[e1,1]]*u[b1,TRIANG[e1,2]]*(EDGE[e1,2,1]*EDGE[e1,3,1] + EDGE[e1,2,2]*EDGE[e1,3,2]) +
                2*u[b1,TRIANG[e1,1]]*u[b1,TRIANG[e1,3]]*(EDGE[e1,2,1]*EDGE[e1,1,1] + EDGE[e1,2,2]*EDGE[e1,1,2]) +
                2*u[b1,TRIANG[e1,2]]*u[b1,TRIANG[e1,3]]*(EDGE[e1,1,1]*EDGE[e1,3,1] + EDGE[e1,1,2]*EDGE[e1,3,2])
            )
        )
    end
    for b1 in 1:dom.BREAK+2
        energy[b1] = sum(integral[b1, e1] for e1 in 1:dom.ELEM)
    end
end

function _initial_position!(problem::PDEProblem, d::PDEDiscretizationDomain, niter)
    ubar = zeros(d.NODES)
    u0 = zeros(d.BREAK+2, d.NODES)
    energy0 = zeros(d.BREAK+2)
    integral0 = zeros(d.BREAK+2, d.ELEM)

    # Calculate an ending point on the other side of the barrier
    for i in 1:niter
        for b1 in 1:d.BREAK+2, n in 1:d.NODES
            u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*d.US[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
        end
        _update_values!(integral0, energy0, u0, problem, d)
        if energy0[end] < 0
            break
        end
        d.UE .*= 2.0
    end

    # Backtrack to the barrier to get a better representation
    ubar .= d.US
    for i in 1:niter
        for b1 in 1:d.BREAK+2, n in 1:d.NODES
            u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*ubar[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
        end
        _update_values!(integral0, energy0, u0, problem, d)
        ebar = maximum(energy0[b1] for b1 in 1:d.BREAK+2 if energy0[b1] < 0.0)
        bbar = minimum(b1 for b1 in 1:d.BREAK+2 if energy0[b1] == ebar)
        ubar .= u0[bbar-1, :]
        d.UE .= u0[bbar, :]
    end

    for b1 in 1:d.BREAK+2, n in 1:d.NODES
        u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*d.US[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
    end
    _update_values!(integral0, energy0, u0, problem, d)
    z0 = maximum(energy0[b1] for b1 in 1:d.BREAK+2)
    return (
        u=u0, energy=energy0, integral=integral0, z=z0,
    )
end

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
function dirichlet_model(nh, backend::AbstractModelerBackend)
    dom = PDEDiscretizationDomain(nh, CIRCLE_DOMAIN)
    pb = PDEProblem(
        0.01,
        fill(1.0, dom.NODES),
        fill(1.0, dom.NODES),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return transition_state_model(pb, dom, backend)
end

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
function henon_model(nh, backend::AbstractModelerBackend)
    dom = PDEDiscretizationDomain(nh, CIRCLE_REC_DOMAIN)
    pb = PDEProblem(
        1.0,
        fill(0.0, dom.NODES),
        sqrt.(dom.COORDS[:, 1].^2 .+ dom.COORDS[:, 2].^2),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return transition_state_model(pb, dom, backend)
end

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
function lane_emden_model(nh, backend::AbstractModelerBackend)
    dom = PDEDiscretizationDomain(nh, RECTANGLE_DOMAIN)
    pb = PDEProblem(
        1.0,
        fill(0.0, dom.NODES),
        fill(1.0, dom.NODES),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return transition_state_model(pb, dom, backend)
end



