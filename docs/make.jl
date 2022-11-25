using SpikingRasters
using Documenter

DocMeta.setdocmeta!(SpikingRasters, :DocTestSetup, :(using SpikingRasters); recursive=true)

makedocs(;
    modules=[SpikingRasters],
    authors="Dylan Festa <dylan.festa@gmail.com>",
    repo="https://github.com/dylanfesta/SpikingRasters.jl/blob/{commit}{path}#{line}",
    sitename="SpikingRasters.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dylanfesta.github.io/SpikingRasters.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dylanfesta/SpikingRasters.jl",
    devbranch="main",
)
