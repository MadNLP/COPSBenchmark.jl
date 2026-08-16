using Documenter
using COPSBenchmark

makedocs(
    sitename = "COPSBenchmark.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/madsuite-org/COPSBenchmark.jl.git",
    devbranch = "main",
    push_preview = true,
)
