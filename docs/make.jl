push!(LOAD_PATH,"../src","../src/wave1D","../src/nbody","../src/TOV","../src/maxwell")

using numex
using Maxwell
using Wave1D
using NBody
using TOV
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
        "Examples" => Any[
            "Euler's Method" => "euler.md",
            "Simple Oscillator" => "sho.md",
            "Nonlinear Oscillator" => "vanderpol.md",
            "N-Body Gravity" => "nbody.md",
            "Stars" => "tov.md",
            "Wave Equation" => "wave1D.md",
            "Maxwell" => "maxwell.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/davidneilsen/numex.jl",
    devbranch="main",
)
