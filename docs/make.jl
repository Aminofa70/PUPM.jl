using PUPM
using Documenter

DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

makedocs(;
    modules=[PUPM],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="PUPM.jl",
    format=Documenter.HTML(;
        canonical="https://Aminofa70.github.io/PUPM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Aminofa70/PUPM.jl",
    devbranch="main",
)
