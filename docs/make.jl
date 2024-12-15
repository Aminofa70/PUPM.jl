
using PUPM
using Documenter

# Set up documentation metadata for doctests
DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

# Build the documentation
makedocs(
    modules = [PUPM],
    authors = "Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo = "https://github.com/Aminofa70/PUPM.jl",
    sitename = "PUPM.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true"
    ),
    pages = [
        "Home" => "index.md",
        # Add more pages as needed
    ]
)

# Deploy the documentation
deploydocs(
    repo = "https://github.com/Aminofa70/PUPM.jl",
    devbranch = "main",
    branch = "gh-pages",
    dirname = "dev",  # Replaces `folder`
    push_preview = false  # Set to true if previews are required
)





# using PUPM
# using Documenter

# DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

# makedocs(;
#     modules=[PUPM],
#     authors= "Aminofa70 <amin.alibakhshi@upm.es> and contributors",
#     repo = "github.com/Aminofa70/PUPM.jl",
#     sitename= "PUPM.jl",
#     format = Documenter.HTML(;
#         prettyurls= get(ENV, "CI", "false") == "true",
#         canonical= "https://Aminofa70.github.io/PUPM.jl",
#         # edit_link= "main",
#         # assets=String[],
#     ),
#     pages=[
#         "Home" => "index.md", 
#         #"Functions" => "functions.md",
#     ],
#    # checkdocs = :none,  # This replaces strict=false# ignores my codes 
# )

# deploydocs(;
#     repo= "github.com/Aminofa70/PUPM.jl",
#     devbranch ="main",
#     branch = "gh-pages",
#     folder = "dev",
#     push_preview = false,
# )
