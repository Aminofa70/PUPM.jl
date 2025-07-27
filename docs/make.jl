# using PUPM
# using Documenter
# using Documenter.Remotes

# # Set up documentation metadata for doctests
# DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

# makedocs(
#     modules = [PUPM],
#     authors = "Aminofa70 <amin.alibakhshi@upm.es> and contributors",
#     repo = Remotes.GitHub("Aminofa70", "PUPM.jl"),
#     sitename = "PUPM.jl",
#     format = Documenter.HTML(
#         prettyurls = get(ENV, "CI", "false") == "true"
#     ),
#     pages = [
#         "Home" => "index.md",
#         "Install" => "install.md"
#     ]
# )

# deploydocs(
#     repo = "https://github.com/Aminofa70/PUPM.jl",
#     branch = "gh-pages",
#     dirname = "",  # or "dev" if you prefer
#     devbranch = "main",
#     push_preview = false
# )


using PUPM
using Documenter
using Documenter.Remotes

DocMeta.setdocmeta!(PUPM, :DocTestSetup, :(using PUPM); recursive=true)

makedocs(
    modules = [PUPM],
    sitename = "PUPM.jl",
    authors = "Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo = Remotes.GitHub("Aminofa70", "PUPM.jl"),  # ✅ Correct usage
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://github.com/Aminofa70/PUPM.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Install" => "install.md",
        "Functions" => "functions.md",
    ]
)

# ✅ Pass devbranch here, NOT in makedocs
deploydocs(
    repo = "github.com/Aminofa70/PUPM.jl",
    devbranch = "main",
    push_preview = false
)
