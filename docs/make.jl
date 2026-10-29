using Lincege
using Documenter

include("example_scripts.jl")
example_scripts = ExampleScripts.generate(joinpath(@__DIR__, "src", "examples"), joinpath(@__DIR__, "..", "examples"))
ExampleScripts.run_examples(example_scripts, @__DIR__)

DocMeta.setdocmeta!(Lincege, :DocTestSetup, :(using Lincege); recursive=true)

makedocs(;
        modules=[Lincege],
        authors="Pranav Seetharaman <pjseetha@uwaterloo.ca> and contributors",
        sitename="LiNCEGE.jl",
        format=Documenter.HTML(;
                canonical="https://rightleftspin.github.io/LiNCEGE.jl",
                edit_link="main",
                assets=String[],
                prettyurls=false,
        ),
        pages=[
                "Home" => "index.md",
                "Examples" => [
                        "Square Lattice Ising Model" => "examples/square_lattice.md",
                        "Square Lattice Symmetric Clusters" => "examples/square_symmetries.md",
                        "Finite Square Lattice (Periodic)" => "examples/finite_square_lattice.md",
                        "Kagome Lattice" => "examples/kagome_lattice.md",
                        "Triangular Lattice" => "examples/triangular_lattice.md",
                        "Lieb Lattice" => "examples/lieb_lattice.md",
                        "One Fifth Depleted Lattice" => "examples/one_fifth_depleted_lattice.md",
                        "Shastry-Sutherland Lattice" => "examples/shastry-sutherland.md",
                        "Pyrochlore Lattice Unit Cell" => "examples/pyrochlore_unit_cell_expansion.md",
                        "Square Lattice Square Expansion" => "examples/square_cluster.md",
                ],
                "Output Formats" => "output_formats.md",
                "Contributing" => "contributing.md",
                "API Reference" => "api.md",
        ],
)

deploydocs(;
        repo="github.com/rightleftspin/LiNCEGE.jl",
        devbranch="main",
)
