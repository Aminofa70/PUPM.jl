using PUPM
using Documenter

DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

makedocs(;
    modules=[PUPM],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo = "github.com/Aminofa70/PUPM.jl",
    sitename= "PUPM.jl",
    format=Documenter.HTML(;
        prettyurls= get(ENV, "CI", "false") == "true",
        canonical= "https://Aminofa70.github.io/PUPM.jl",
        # edit_link= "main",
        # assets=String[],
    ),
    pages=[
        "Home" => "index.md", 
        "Functions" => "functions.md",
    ],
   # checkdocs = :none,  # This replaces strict=false# ignores my codes 
)

deploydocs(;
    repo="github.com/Aminofa70/PUPM.jl",
    devbranch="main",
    push_preview = false,
)
