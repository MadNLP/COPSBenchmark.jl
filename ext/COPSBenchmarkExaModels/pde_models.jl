
@inline function COPSBenchmark.transition_state_model(::ExaModelsBackend, problem, dom::COPSBenchmark.PDEDiscretizationDomain; T = Float64, backend = nothing, kwargs...)
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    # Precompute T-typed constants so the kernel expressions stay fp32-clean:
    # the scalar `a` and the integer divisors `/2`, `1/(DIMEN+1)` would
    # otherwise inject Float64 into the GPU AD kernels.
    aT = T(a)
    half = T(1) / T(2)
    inv_dim = T(1) / T(dom.DIMEN + 1)
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
    core = ExaModels.ExaCore(T; backend = backend, concrete = Val(true))

    ExaModels.@add_var(core, u, 1:dom.BREAK+2, 1:dom.NODES; start=x0.u)
    ExaModels.@add_var(core, integral, 1:dom.BREAK+2, 1:dom.ELEM)
    ExaModels.@add_var(core, z, 1; start=x0.z)

    ExaModels.@add_obj(core, z[1])

    ExaModels.@add_con(
        core,
        c1,
        - z[1] for b1 in 1:dom.BREAK + 2;
        lcon = -Inf,
        ucon = 0.0,
    )

    ExaModels.@add_con!(
        core,
        c1,
        b1 => integral[b1, e1] for b1 in 1:dom.BREAK +2, e1 in 1:dom.ELEM
    )

    ExaModels.@add_con(
        core,
        c2,
        dom.BREAK+1;
        lcon=-Inf,
        ucon=H^2,
    )

    ExaModels.@add_con!(
        core,
        c2,
        b1 => (u[b1+1, n] - u[b1, n])^2 for b1 in 1:dom.BREAK+1, n in 1:dom.NODES
    )

    ExaModels.@add_con(
        core,
        c3,
        AREA*(
            aT / (8*AREA^2)*(
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

    ExaModels.@add_con!(
        core,
        c3,
        (b1, e1) => AREA * inv_dim *
                    (b*u[b1,TRIANG]^2*half - c*u[b1,TRIANG]^(p+1)/(p+1)+ d*u[b1, TRIANG])
                    for b1 in 1:dom.BREAK+2, (e1, AREA, b, c, d, p, TRIANG) in array2
    )

    # Boundary
    boundary_nodes = findall(isequal(1), dom.BNDRY)
    ExaModels.@add_con(
        core,
        c4,
        u[b1+1, n] for b1 in 1:dom.BREAK, n in boundary_nodes
    )
    ExaModels.@add_con(
        core,
        c5,
        u[1, n] for n in 1:dom.NODES;
        lcon=dom.US,
        ucon=dom.US,
    )
    ExaModels.@add_con(
        core,
        c6,
        u[dom.BREAK+2, n] for n in 1:dom.NODES;
        lcon=dom.UE,
        ucon=dom.UE,
    )
    return ExaModels.ExaModel(core; kwargs...)
end


# transition_state family: one deterministic init from the CURRENT domain
# state. The eager `_initial_position!` MUTATES `dom.UE` in place — and `dom`
# aliases the global `*_DOMAIN` constants — so each init here runs on a
# deepcopy: the lazy model computes exactly what the eager constructor (run
# next by the harness, from the same global state) computes, without advancing
# that state itself. `UE` (post-run) is returned because the eager captures it
# for the c6 bounds and `H` AFTER the mutation (see transition_state_model above).
function _pde_lazy_init(rawdom, mk_pb, nh)
    dom = COPSBenchmark.PDEDiscretizationDomain(nh, deepcopy(rawdom))
    pb = mk_pb(dom)
    x0 = COPSBenchmark._initial_position!(pb, dom, 10)
    ALPHA = 2.0
    H = ALPHA / (dom.BREAK + 1) * sqrt(sum((dom.US[n] - dom.UE[n])^2 for n in 1:dom.NODES))
    return (u0 = x0.u, z0 = x0.z, H2 = H^2, UE = dom.UE)
end

# The three problem instances (COPS/src/pde_utils.jl:159-215); the PDEProblem
# coefficient vectors depend only on the mesh (NODES/COORDS), never on nh.
_pde_pb_dirichlet(dom) = COPSBenchmark.PDEProblem(
    0.01, fill(1.0, dom.NODES), fill(1.0, dom.NODES), fill(0.0, dom.NODES), fill(3.0, dom.NODES))
_pde_pb_henon(dom) = COPSBenchmark.PDEProblem(
    1.0, fill(0.0, dom.NODES), sqrt.(dom.COORDS[:, 1] .^ 2 .+ dom.COORDS[:, 2] .^ 2),
    fill(0.0, dom.NODES), fill(3.0, dom.NODES))
_pde_pb_lane_emden(dom) = COPSBenchmark.PDEProblem(
    1.0, fill(0.0, dom.NODES), fill(1.0, dom.NODES), fill(0.0, dom.NODES), fill(3.0, dom.NODES))

# 1:1 transcription of transition_state_model (pde_models.jl). `nh` is either
# `args.nh` (dirichlet/henon/lane_emden) or a literal Int (transition_state);
# u0/z0/H2/UE fill the corresponding value slots (Deferred or concrete).
# array1/array2/AREA/EDGE/BNDRY/US depend only on the mesh, so they are
# computed concretely here (BREAK is unused; nh = 0 placeholder).
function _lazy_transition_state(nh, rawdom, mk_pb, u0, z0, H2, UE)
    dom = COPSBenchmark.PDEDiscretizationDomain(0, rawdom)
    pb = mk_pb(dom)
    b, c, d, p = pb.b, pb.c, pb.d, pb.p
    aT = pb.a
    half = 0.5
    inv_dim = 1.0 / (dom.DIMEN + 1)
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

    core = ExaCore(concrete = Val(true))

    core, u = add_var(core, 1:nh+2, 1:dom.NODES; start = u0)
    core, integral = add_var(core, 1:nh+2, 1:dom.ELEM)
    core, z = add_var(core, 1; start = z0)

    core, _ = add_obj(core, z[1] for _ in 1:1)

    core, c1 = add_con(core, (-z[1] for b1 in 1:nh+2); lcon = -Inf, ucon = 0.0)
    core, _ = add_con!(core, c1, b1 => integral[b1, e1] for b1 in 1:nh+2, e1 in 1:dom.ELEM)

    # empty-dims constraint: spelled 1:nh+1 — a bare ArgLike dim would
    # dispatch to add_con's scalar-expression form instead.
    core, c2 = add_con(core, 1:nh+1; lcon = -Inf, ucon = H2)
    core, _ = add_con!(core, c2, b1 => (u[b1+1, n] - u[b1, n])^2 for b1 in 1:nh+1, n in 1:dom.NODES)

    core, c3 = add_con(core,
        AREA*(
            aT / (8*AREA^2)*(
                u[b1,TRIANG1]^2*(EDGE_21^2 + EDGE_22^2) +
                u[b1,TRIANG2]^2*(EDGE_31^2 + EDGE_32^2) +
                u[b1,TRIANG3]^2*(EDGE_11^2 + EDGE_12^2) +
                2*u[b1,TRIANG1]*u[b1,TRIANG2]*(EDGE_21*EDGE_31 + EDGE_22*EDGE_32) +
                2*u[b1,TRIANG1]*u[b1,TRIANG3]*(EDGE_21*EDGE_11 + EDGE_22*EDGE_12) +
                2*u[b1,TRIANG2]*u[b1,TRIANG3]*(EDGE_11*EDGE_31 + EDGE_12*EDGE_32)
            )
        )
        - integral[b1, e1]
        for b1 in 1:nh+2, (e1, TRIANG1, TRIANG2, TRIANG3, AREA, EDGE_11, EDGE_12, EDGE_21, EDGE_22, EDGE_31, EDGE_32) in array1)
    core, _ = add_con!(core, c3,
        (b1, e1) => AREA * inv_dim *
                    (b*u[b1,TRIANG]^2*half - c*u[b1,TRIANG]^(p+1)/(p+1) + d*u[b1, TRIANG])
        for b1 in 1:nh+2, (e1, AREA, b, c, d, p, TRIANG) in array2)

    # Boundary
    boundary_nodes = findall(isequal(1), dom.BNDRY)
    core, _ = add_con(core, u[b1+1, n] for b1 in 1:nh, n in boundary_nodes)
    core, _ = add_con(core, (u[1, n] for n in 1:dom.NODES); lcon = dom.US, ucon = dom.US)
    core, _ = add_con(core, (u[nh+2, n] for n in 1:dom.NODES); lcon = UE, ucon = UE)
    core
end

# transition_state-family builder over a deferred nh: every UE-derived value
# slot is a Deferred running `_pde_lazy_init` at the resolved nh.
function _lazy_pde_core(nh, rawdom, mk_pb)
    return _lazy_transition_state(nh, rawdom, mk_pb,
        Deferred(a -> _pde_lazy_init(rawdom, mk_pb, a.nh).u0),
        Deferred(a -> _pde_lazy_init(rawdom, mk_pb, a.nh).z0),
        Deferred(a -> _pde_lazy_init(rawdom, mk_pb, a.nh).H2),
        Deferred(a -> _pde_lazy_init(rawdom, mk_pb, a.nh).UE))
end

function COPSBenchmark.dirichlet_core()
    args = ExaModels.ArgTracer()
    return _lazy_pde_core(args.nh, COPSBenchmark.CIRCLE_DOMAIN, _pde_pb_dirichlet)
end

function COPSBenchmark.henon_core()
    args = ExaModels.ArgTracer()
    return _lazy_pde_core(args.nh, COPSBenchmark.CIRCLE_REC_DOMAIN, _pde_pb_henon)
end

function COPSBenchmark.lane_emden_core()
    args = ExaModels.ArgTracer()
    return _lazy_pde_core(args.nh, COPSBenchmark.RECTANGLE_DOMAIN, _pde_pb_lane_emden)
end

# Fixed-size instance: the canonical dirichlet problem on CIRCLE_DOMAIN at
# nh = 5 (the eager constructor takes prebuilt problem/domain structures, not
# a size). Value slots are computed eagerly at build time, from a deepcopy.
function COPSBenchmark.transition_state_core()
    init = _pde_lazy_init(COPSBenchmark.CIRCLE_DOMAIN, _pde_pb_dirichlet, 5)
    return _lazy_transition_state(5, COPSBenchmark.CIRCLE_DOMAIN, _pde_pb_dirichlet,
        init.u0, init.z0, init.H2, init.UE)
end
