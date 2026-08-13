module COPSBenchmarkExaModels

using Random
import COPSBenchmark
import COPSBenchmark: ExaModelsBackend
using ExaModels

# ── Deferred callables in recipes ─────────────────────────────────────────────
#
# Every function a recipe defers — start generators, index-set builders, data
# tables — must be NAMED and owned by THIS PACKAGE: not an anonymous closure,
# not an extension-owned function, not anything in `Main`.  Two independent
# mechanisms require it.  The serialized core carries the function's TYPE,
# which another process can resolve only through the owning package's `PkgId`;
# and an AOT-compiled library calls argument functions by NAME, which only a
# package import makes reachable.  `Base.Fix2` over such a function is fine —
# the wrapper's type parameters stay named.

include("bearing.jl")
include("camshape.jl")
include("catmix.jl")
include("chain.jl")
include("elec.jl")
include("gasoil.jl")
include("glider.jl")
include("marine.jl")
include("methanol.jl")
include("minsurf.jl")
include("pinene.jl")
include("robot.jl")
include("rocket.jl")
include("steering.jl")
include("channel.jl")
include("torsion.jl")
include("polygon.jl")
include("triangle.jl")
include("tetra.jl")
include("pde_models.jl")

end

