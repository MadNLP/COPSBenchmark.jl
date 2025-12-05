
function COPSBenchmark.transition_state_model(problem, dom::COPSBenchmark.PDEDiscretizationDomain, ::ExaModelsBackend; T = Float64, backend = nothing, kwargs...)
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    x0 = COPSBenchmark._initial_position!(problem, dom, 10)
    array1 = [
        (
            e1,
            dom.TRIANG[e1, 1],
            dom.TRIANG[e1, 2],
            dom.TRIANG[e1, 3],
            dom.AREA[e1],
            dom.EDGE[e1, 1, 1],
            dom.EDGE[e1, 1, 2],
            dom.EDGE[e1, 2, 1],
            dom.EDGE[e1, 2, 2],
            dom.EDGE[e1, 3, 1],
            dom.EDGE[e1, 3, 2]
        ) for e1 in 1:dom.ELEM
    ]
    array2 = [
    (
        e1 = e1,
        AREA = dom.AREA[e1],
        b = b[dom.TRIANG[e1, c1]],
        c = c[dom.TRIANG[e1, c1]],
        d = d[dom.TRIANG[e1, c1]],
        p = p[dom.TRIANG[e1, c1]],
        TRIANG = dom.TRIANG[e1, c1]
    )
    for e1 in 1:dom.ELEM for c1 in 1:(dom.DIMEN + 1)
    ]
    #_proto_model()

    ALPHA = 2.0
    H = ALPHA / (dom.BREAK+1) * sqrt(sum((dom.US[n] - dom.UE[n])^2 for n in 1:dom.NODES))

    # Build optimization problem
    core = ExaModels.ExaCore(T; backend = backend)

    u = ExaModels.variable(core, 1:dom.BREAK+2, 1:dom.NODES; start=x0.u)
    integral = ExaModels.variable(core, 1:dom.BREAK+2, 1:dom.ELEM)
    z = ExaModels.variable(core, 1; start=x0.z)

    ExaModels.objective(core, z[1])

    c1 = ExaModels.constraint(
        core,
        - z[1] for b1 in 1:dom.BREAK + 2;
        lcon = -Inf,
        ucon = 0.0,
    )

    ExaModels.constraint!(
        core,
        c1,
        b1 => integral[b1, e1] for b1 in 1:dom.BREAK +2, e1 in 1:dom.ELEM
    )

    c2 = ExaModels.constraint(
        core,
        dom.BREAK+1;
        lcon=-Inf,
        ucon=H^2,
    )

    ExaModels.constraint!(
        core,
        c2,
        b1 => (u[b1+1, n] - u[b1, n])^2 for b1 in 1:dom.BREAK+1, n in 1:dom.NODES
    )

    c3 = ExaModels.constraint(
        core,
        AREA*(
            a / (8*AREA^2)*(
                u[b1,TRIANG1]^2*(EDGE_21^2 + EDGE_22^2) +
                u[b1,TRIANG2]^2*(EDGE_31^2 + EDGE_32^2) +
                u[b1,TRIANG3]^2*(EDGE_11^2 + EDGE_12^2) +
                2*u[b1,TRIANG1]*u[b1,TRIANG2]*(EDGE_21*EDGE_31 + EDGE_22*EDGE_32) +
                2*u[b1,TRIANG1]*u[b1,TRIANG3]*(EDGE_21*EDGE_11 + EDGE_22*EDGE_12) +
                2*u[b1,TRIANG2]*u[b1,TRIANG3]*(EDGE_11*EDGE_31 + EDGE_12*EDGE_32)
                )
            )
            - integral[b1, e1]
        for b1 in 1:dom.BREAK+2, (e1, TRIANG1, TRIANG2, TRIANG3, AREA, EDGE_11, EDGE_12, EDGE_21, EDGE_22, EDGE_31, EDGE_32) in array1
    )

    ExaModels.constraint!(
        core,
        c3,
        (b1, e1) => AREA* 1 / (dom.DIMEN+1) *
                    (b*u[b1,TRIANG]^2/2- c*u[b1,TRIANG]^(p+1)/(p+1)+ d*u[b1, TRIANG])
                    for b1 in 1:dom.BREAK+2, (e1, AREA, b, c, d, p, TRIANG) in array2
    )

    # Boundary
    boundary_nodes = findall(isequal(1), dom.BNDRY)
    ExaModels.constraint(
        core,
        u[b1+1, n] for b1 in 1:dom.BREAK, n in boundary_nodes
    )
    ExaModels.constraint(
        core,
        u[1, n] for n in 1:dom.NODES;
        lcon=dom.US,
        ucon=dom.US,
    )
    ExaModels.constraint(
        core,
        u[dom.BREAK+2, n] for n in 1:dom.NODES;
        lcon=dom.UE,
        ucon=dom.UE,
    )
    return ExaModels.ExaModel(core; kwargs...)
end

