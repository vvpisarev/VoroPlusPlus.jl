using Documenter
using VoroPlusPlus

makedocs(;
    modules = [VoroPlusPlus],
    checkdocs = :exports,
    sitename = "VoroPlusPlus.jl documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    repo = Remotes.GitHub("vvpisarev", "VoroPlusPlus.jl"),
    pages = [
        "Home" => "index.md",
        "Containers" => "container.md",
        "Voronoi cells" => "cell.md",
        "Walls" => "walls.md",
        "Iteration over container" => "iteration.md",
        # "Release Notes" => "release_notes.md"
    ],
    warnonly = true,
)

deploydocs(;
  repo = "github.com/vvpisarev/VoroPlusPlus.jl",
  devurl = "dev",
  devbranch = "main",
  versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
