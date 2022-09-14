using numex
using Documenter

DocMeta.setdocmeta!(numex, :DocTestSetup, :(using numex); recursive=true)

makedocs(;
    modules=[numex],
    authors="David Neilsen <david.neilsen@byu.edu> and contributors",
    repo="https://github.com/davidneilsen/numex.jl/blob/{commit}{path}#{line}",
    sitename="numex.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://davidneilsen.github.io/numex.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/davidneilsen/numex.jl",
    devbranch="main",
)
