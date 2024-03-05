module COPSBenchmark

using Random
using JuMP


include("bearing.jl")
include("chain.jl")
include("camshape.jl")
include("catmix.jl")
include("channel.jl")
include("elec.jl")
include("gasoil.jl")
include("glider.jl")
include("marine.jl")
include("methanol.jl")
include("pinene.jl")
include("polygon.jl")
include("robot.jl")
include("rocket.jl")
include("steering.jl")
include("tetra.jl")
include("torsion.jl")
include("triangle.jl")

include("pde_utils.jl")
include("dirichlet.jl")
include("henon.jl")
include("lane_emden.jl")

end # module COPSBenchmark
