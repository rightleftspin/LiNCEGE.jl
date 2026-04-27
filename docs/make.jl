using LiNCEGE
using Documenter

DocMeta.setdocmeta!(LiNCEGE, :DocTestSetup, :(using LiNCEGE); recursive=true)

makedocs(;
        modules=[LiNCEGE],
        authors="Pranav Seetharaman <pranav@myrdd.info> and contributors",
        sitename="LiNCEGE.jl",
        format=Documenter.HTML(;
                canonical="https://rightleftspin.github.io/LiNCEGE.jl",
                edit_link="main",
                assets=String[],
        ),
        pages=[
                "Home" => "index.md",
                "Examples" => [
                        "Square Lattice" => "examples/square_lattice.md",
                        "Kagome Lattice" => "examples/kagome_lattice.md",
                        "Pyrochlore Lattice Unit Cell" => "examples/pyrochlore_unit_cell_expansion.md",
                        "Square Lattice Cluster" => "examples/square_cluster.md",
                ],
        ],
)

deploydocs(;
        repo="github.com/rightleftspin/LiNCEGE.jl",
        devbranch="main",
)
