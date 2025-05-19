using BenchmarkTools
using NLCE 

begin # SETUP
    square_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[1, 0], [0, 1]])
    triangular_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[0.5, sqrt(3) / 2], [1, 0]])
    kagome_lattice = Dict(
            "Basis" => [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]],
            "Primitive Vectors" => [[2, 0], [1, sqrt(3)]],
        )

    const SUITE = BenchmarkGroup()
    nlce = SUITE["NLCE"] = BenchmarkGroup()
end

begin # Lattice Creation
    nlce["creation"] = BenchmarkGroup()

    neighbors = [1]
    max_order = 6
    
    nlce["creation"]["create square lattice"] = @benchmarkable NLCE.SiteExpansionBundle(
            square_lattice["Basis"],
            square_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )
    nlce["creation"]["create triangular lattice"] = @benchmarkable NLCE.SiteExpansionBundle(
            triangular_lattice["Basis"],
            triangular_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )
    nlce["creation"]["create kagome lattice"] = @benchmarkable NLCE.SiteExpansionBundle(
            kagome_lattice["Basis"],
            kagome_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )
end

begin # Evaluation
    res = BenchmarkTools.run(SUITE, verbose=true)
    median_nlce = median(res["NLCE"])
end

