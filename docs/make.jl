using WignerSeitz
using Documenter

DocMeta.setdocmeta!(WignerSeitz, :DocTestSetup, :(using WignerSeitz); recursive=true)

makedocs(;
    modules=[WignerSeitz],
    authors="Brandon Flores",
    repo="https://github.com/brainandforce/WignerSeitz.jl/blob/{commit}{path}#{line}",
    sitename="WignerSeitz.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brainandforce.github.io/WignerSeitz.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brainandforce/WignerSeitz.jl",
    devbranch="main",
)
