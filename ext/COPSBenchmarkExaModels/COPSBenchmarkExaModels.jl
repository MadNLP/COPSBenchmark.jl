module COPSBenchmarkExaModels

using Random
import COPSBenchmark
import COPSBenchmark: ExaModelsBackend
using ExaModels

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
include("pde_models.jl")

end

