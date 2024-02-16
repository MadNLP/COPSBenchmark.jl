module COPSBenchmark

using Random
using JuMP


include("bearing.jl")
include("camshape.jl")
include("elec.jl")
include("gasoil.jl")
include("marine.jl")
include("pinene.jl")
include("robot.jl")
include("rocket.jl")
include("steering.jl")

include("pde_utils.jl")
include("dirichlet.jl")
include("henon.jl")
include("lane_emden.jl")

end # module COPSBenchmark
