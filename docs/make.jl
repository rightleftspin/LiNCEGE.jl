using LINCEGE
using Documenter

DocMeta.setdocmeta!(LINCEGE, :DocTestSetup, :(using LINCEGE); recursive=true)

makedocs(;
        modules=[LINCEGE],
        authors="Pranav Seetharaman <pranav@myrdd.info> and contributors",
        sitename="LINCEGE.jl",
        format=Documenter.HTML(;
                canonical="https://rightleftspin.github.io/LINCEGE.jl",
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
        repo="github.com/rightleftspin/LINCEGE.jl",
        devbranch="main",
)
