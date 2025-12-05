
function COPSBenchmark.transition_state_model(problem, dom::COPSBenchmark.PDEDiscretizationDomain, ::JuMPBackend)
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    x0 = COPSBenchmark._initial_position!(problem, dom, 10)

    ALPHA = 2.0
    H = ALPHA / (dom.BREAK+1) * sqrt(sum((dom.US[n] - dom.UE[n])^2 for n in 1:dom.NODES))

    model = Model()

    @variable(model, u[b1=0:dom.BREAK+1, n=1:dom.NODES], start=x0.u[b1+1, n])
    @variable(model, z, start=x0.z)

    @objective(model, Min, z)

    @expression(
        model,
        integral[b1 in 0:dom.BREAK+1, e1 in 1:dom.ELEM],
        dom.AREA[e1]*(
        1 / (dom.DIMEN+1) *
            (sum((b[dom.TRIANG[e1,c1]]*u[b1,dom.TRIANG[e1,c1]]^2/2-
                        c[dom.TRIANG[e1,c1]]*u[b1,dom.TRIANG[e1,c1]]^(p[dom.TRIANG[e1,c1]]+1)/(p[dom.TRIANG[e1,c1]]+1)+
                        d[dom.TRIANG[e1,c1]]*u[b1,dom.TRIANG[e1,c1]]) for c1 in 1:dom.DIMEN+1)) +
        a / (8*dom.AREA[e1]^2)*(
            u[b1,dom.TRIANG[e1,1]]^2*(dom.EDGE[e1,2,1]^2 + dom.EDGE[e1,2,2]^2) +
            u[b1,dom.TRIANG[e1,2]]^2*(dom.EDGE[e1,3,1]^2 + dom.EDGE[e1,3,2]^2) +
            u[b1,dom.TRIANG[e1,3]]^2*(dom.EDGE[e1,1,1]^2 + dom.EDGE[e1,1,2]^2) +
            2*u[b1,dom.TRIANG[e1,1]]*u[b1,dom.TRIANG[e1,2]]*(dom.EDGE[e1,2,1]*dom.EDGE[e1,3,1] + dom.EDGE[e1,2,2]*dom.EDGE[e1,3,2]) +
            2*u[b1,dom.TRIANG[e1,1]]*u[b1,dom.TRIANG[e1,3]]*(dom.EDGE[e1,2,1]*dom.EDGE[e1,1,1] + dom.EDGE[e1,2,2]*dom.EDGE[e1,1,2]) +
            2*u[b1,dom.TRIANG[e1,2]]*u[b1,dom.TRIANG[e1,3]]*(dom.EDGE[e1,1,1]*dom.EDGE[e1,3,1] + dom.EDGE[e1,1,2]*dom.EDGE[e1,3,2])
            )
        )
    )
    @expression(model, energy[b1 in 0:dom.BREAK+1], sum(integral[b1, e1] for e1 in 1:dom.ELEM))

    @constraint(model, [b1 in 1:dom.BREAK], z >= energy[b1])
    @constraint(model, [b1 in 0:dom.BREAK], sum((u[b1+1, n] - u[b1, n])^2 for n in 1:dom.NODES) <= H^2)

    # Boundary
    boundary_nodes = findall(isequal(1), dom.BNDRY)
    @constraint(
        model,
        [b1 in 1:dom.BREAK, n in boundary_nodes],
        u[b1, n] == 0.0,
    )
    @constraint(model, [n in 1:dom.NODES], u[0, n] == dom.US[n])
    @constraint(model, [n in 1:dom.NODES], u[dom.BREAK+1, n] == dom.UE[n])

    return model
end

