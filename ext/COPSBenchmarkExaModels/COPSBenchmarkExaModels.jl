module COPSBenchmarkExaModels

using Random
import COPSBenchmark
import COPSBenchmark: ExaModelsBackend
using ExaModels

# `Deferred(f)` (not exported): a function of the resolved args, evaluated
# once per materialization — used by the *_core builders as a value slot
# (computed start arrays, bounds) or as an index set (pars records carrying
# precomputed index data per element).
const Deferred = ExaModels.Deferred

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

