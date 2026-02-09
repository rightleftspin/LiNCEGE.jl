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
    ],
)

deploydocs(;
    repo="github.com/rightleftspin/LINCEGE.jl",
    devbranch="main",
)
