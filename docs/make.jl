using Documenter
using COPSBenchmark

makedocs(
    sitename = "COPSBenchmark.jl",
    modules = [COPSBenchmark],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/MadNLP/COPSBenchmark.jl.git",
    devbranch = "main",
    push_preview = true,
)
