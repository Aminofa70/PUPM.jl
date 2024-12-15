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
    repo = "https://github.com/Aminofa70/PUPM.jl",  # Full HTTPS URL to your repository
    branch = "gh-pages",                           # GitHub Pages branch
   # dirname = "stable",                            # Folder for stable docs
   dirname = "",                                # Folder for development docs
    devbranch = "main",                            # Development branch
    push_preview = false                           # Optional: Disable previews
)





