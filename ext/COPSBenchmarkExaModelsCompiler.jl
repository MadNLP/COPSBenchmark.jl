# Compiling every COPS model into one shared library.
#
# The instantiation spec belongs to the model, not to the call: three of these
# problems are posed on a grid and take two sizes where the other fourteen take
# one, and the collocation problems carry their element type into the deferred
# data functions.  So the caller says how big and in what precision, and this
# extension knows which model wants that as `(n,)` and which as `(nx, ny)`.
#
# It lives in an extension so the package never depends on the compiler: a
# benchmark set is a modelling package, and a compiler toolchain is not
# something you should acquire by loading one.
module COPSBenchmarkExaModelsCompiler

import COPSBenchmark as CB
import ExaModelsCompiler
using ExaModelsCompiler: compile_library, select

# Posed on a 2-D grid; everything else takes a single discretization size.
const GRID_MODELS = (:bearing, :minsurf, :torsion)

_problems() = sort!([
    Symbol(chopsuffix(string(n), "_recipe")) for n in Base.names(CB; all = true)
    if endswith(string(n), "_recipe") && !startswith(string(n), "#")
])

function ExaModelsCompiler.compile_all(
    ::Val{CB};
    path = "@cops",
    sizes = 100,
    grid_sizes = (10, 10),
    T = Float64,
    only = nothing,
    exclude = (),
    kwargs...,
)
    b = CB.ExaModelsBackend()
    # Assemble every model, then let `select` filter: it is the shared contract
    # and refuses a name this package does not provide, rather than quietly
    # compiling a library missing the model the caller asked for.  Building a
    # recipe is cheap — the expense is the compile, which only the survivors
    # reach.
    models = map(_problems()) do name
        recipe = getfield(CB, Symbol(name, :_recipe))
        args = getfield(CB, Symbol(name, :_args))
        spec = name in GRID_MODELS ? args(b, grid_sizes...) : args(b, sizes)
        name => (recipe(b; T = T), spec...)
    end
    return compile_library(path, select(models; only, exclude)...; kwargs...)
end

end # module
