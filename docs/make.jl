using Documenter, MendelEstimateFrequencies

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "MendelEstimateFrequencies",
    modules = [MendelEstimateFrequencies],
    authors = "Jeanette Papp",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/MendelEstimateFrequencies.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
