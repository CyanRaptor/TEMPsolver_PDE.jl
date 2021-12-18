using TEMPsolver_PDE
using Documenter

DocMeta.setdocmeta!(TEMPsolver_PDE, :DocTestSetup, :(using TEMPsolver_PDE); recursive=true)

makedocs(;
    modules=[TEMPsolver_PDE],
    authors="Amir Farzin",
    repo="https://github.com/CyanRaptor/TEMPsolver_PDE.jl/blob/{commit}{path}#{line}",
    sitename="TEMPsolver_PDE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
