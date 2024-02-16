module COPSBenchmark

using Random
using JuMP

include("bearing.jl")
include("camshape.jl")
include("catmix.jl")
include("elec.jl")
include("gasoil.jl")
include("marine.jl")
include("pinene.jl")
include("robot.jl")
include("rocket.jl")
include("steering.jl")
include("torsion.jl")

end # module COPSBenchmark
