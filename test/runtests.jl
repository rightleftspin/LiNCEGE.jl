using Test
using Plots
using JSON3
using LINCEGE
import LINCEGE: basis_size, dimension, shift_unit_cell,
        centers, max_order, n_unique_sites, neighbors, bond_matrix

# ---------------------------------------------------------------------------
# Shared unit cell definitions reused across multiple testsets
# ---------------------------------------------------------------------------

square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_symmetries::Vector{Matrix{Float64}} = [[1 0; 0 1], [0 1; -1 0], [-1 0; 0 -1], [0 -1; 1 0], [0 1; 1 0], [0 -1; -1 0], [1 0; 0 -1], [-1 0; 0 1]]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])

kagome_basis = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
kagome_pvecs = [[2.0, 0.0], [1.0, sqrt(3)]]
kagome_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
kagome_uc = UnitCell(kagome_basis, kagome_pvecs, kagome_bonds, [1, 1, 1])

ss_pvecs = [[(sqrt(3) + 1) / sqrt(2), 0.0], [0.0, (sqrt(3) + 1) / sqrt(2)]]
ss_basis = [[0.0, 0.0],
        [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))],
        [sqrt(2) / 2, sqrt(6) / 2],
        [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
ss_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
        Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
        Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)]
ss_uc = UnitCell(ss_basis, ss_pvecs, ss_bonds, [1, 2, 3, 4])

one_fifth_pvecs = [[2.0, 1.0], [-1.0, 2.0]]
one_fifth_basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0], [1.0, 1.0]]
one_fifth_bonds = [Bond(1, 2, [0, 0], 1), Bond(1, 3, [0, -1], 1),
        Bond(2, 4, [0, 0], 1), Bond(3, 4, [0, 1], 1),
        Bond(2, 3, [0, 0], 2), Bond(1, 4, [-1, 0], 2)]
one_fifth_uc = UnitCell(one_fifth_basis, one_fifth_pvecs, one_fifth_bonds, [1, 1, 1, 1])

cube_pvecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
cube_bonds = [Bond(1, 1, [1, 0, 0], 1), Bond(1, 1, [0, 1, 0], 1), Bond(1, 1, [0, 0, 1], 1)]
cube_uc = UnitCell([[0.0, 0.0, 0.0]], cube_pvecs, cube_bonds, [1])


# ---------------------------------------------------------------------------

@testset verbose = true "LINCEGE.jl" begin

        # -------------------------------------------------------------------
        @testset verbose = true "Vertices" begin

                @testset "LatticeVertices" begin
                        lv1 = LatticeVertices([1, 3, 5])
                        lv2 = LatticeVertices([3, 4, 5])

                        @test collect(lv1) == [1, 3, 5]
                        @test sort(lv1) == lv1
                        @test intersect(lv1, lv2) == LatticeVertices([3, 5])
                        @test setdiff(lv1, lv2) == LatticeVertices(1)
                        @test union(lv1, lv2) == LatticeVertices([1, 3, 4, 5])
                        @test (3 in lv1) == true
                        @test eltype(lv1) == Int
                end

                @testset "ExpansionVertices" begin
                        ev1 = ExpansionVertices([2, 4, 6])
                        ev2 = ExpansionVertices([4, 5, 6])

                        @test collect(ev1) == [2, 4, 6]
                        @test sort(ev1) == ev1
                        @test intersect(ev1, ev2) == ExpansionVertices([4, 6])
                        @test setdiff(ev1, ev2) == ExpansionVertices(2)
                        @test union(ev1, ev2) == ExpansionVertices([2, 4, 5, 6])
                        @test (4 in ev1) == true
                        @test eltype(ev1) == Int
                end

                @testset "AbstractVertices" begin
                        lv = LatticeVertices([1, 2])
                        ev = ExpansionVertices([3, 4])

                        @test length(lv) == 2
                        @test contains(lv, 1) == true
                        @test haskey(lv, 1) == true

                        v = [10, 20, 30, 40]
                        m = [10 20 30 40; 50 60 70 80; 90 100 110 120; 130 140 150 160]

                        @test v[lv] == [10, 20]
                        @test m[ev, ev] == [110 120; 150 160]

                        iterate_lv = collect(lv)
                        for (index, value) in enumerate(lv)
                                @test iterate_lv[index] == value
                        end

                        show_output = IOBuffer()
                        show(show_output, lv)
                        @test String(take!(show_output)) == "Vertices: [1, 2]"
                end

        end # Vertices

        # -------------------------------------------------------------------
        @testset verbose = true "UnitCells" begin

                @testset "Site Expansion Structure" begin
                        @test basis_size(square_uc) == 1
                        @test dimension(square_uc) == 2
                        @test shift_unit_cell(square_uc, [1, 0, 1]) == [1.0, 0.0]
                        @test shift_unit_cell(square_uc, [2, -2, 1]) == [2.0, -2.0]
                        @test shift_unit_cell(square_uc, [1 2 3; 0 -2 1; 1 1 1]) == [1.0 2.0 3.0; 0.0 -2.0 1.0]
                end

                @testset "Cluster Expansion Structure" begin end

                @testset "Visualization" begin
                        @test image_unit_cell(square_uc) isa Plots.Plot{Plots.GRBackend}
                        @test image_unit_cell(kagome_uc) isa Plots.Plot{Plots.GRBackend}
                        @test image_unit_cell(one_fifth_uc) isa Plots.Plot{Plots.GRBackend}
                        @test image_unit_cell(ss_uc) isa Plots.Plot{Plots.GRBackend}
                        @test image_unit_cell(cube_uc) isa Plots.Plot{Plots.GRBackend}
                end

        end # UnitCells

        # -------------------------------------------------------------------
        @testset verbose = true "Lattices" begin

                @testset verbose = true "SiteExpansionLattice" begin

                        @testset "Square Lattice" begin
                                lattice = SiteExpansionLattice(2, square_uc)

                                @test centers(lattice) == ExpansionVertices([13])
                                @test max_order(lattice) == UInt8(2)
                                @test n_unique_sites(lattice) == 1
                                @test neighbors(lattice, centers(lattice)) == ExpansionVertices([12, 14, 8, 18])
                        end

                        @testset "Kagome Lattice" begin
                                lattice = SiteExpansionLattice(2, kagome_uc)

                                @test centers(lattice) == ExpansionVertices([13, 38, 63])
                                @test max_order(lattice) == UInt8(2)
                                @test n_unique_sites(lattice) == 3
                                @test neighbors(lattice, LatticeVertices([38])) == LatticeVertices([13, 59, 14, 63])
                        end

                end

                @testset "Bond Type Correctness" begin
                        sq_lat = SiteExpansionLattice(2, square_uc)
                        @test sort(unique(filter(!=(0), bond_matrix(sq_lat)))) == [1]

                        ss_lat = SiteExpansionLattice(2, ss_uc)
                        @test sort(unique(filter(!=(0), bond_matrix(ss_lat)))) == [1, 2, 3]
                end

        end # Lattices

        # -------------------------------------------------------------------
        @testset verbose = true "Clusters" begin

                @testset "Translation clustering (Square)" begin
                        lattice = SiteExpansionLattice(2, square_uc)
                        clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(clusters, lattice)

                        @test length(filter(c -> length(c) == 1, collect(clusters))) == 1
                        @test length(filter(c -> length(c) == 2, collect(clusters))) == 2
                end

                @testset "Symmetric clustering (Square)" begin
                        lattice = SiteExpansionLattice(8, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)

                        sym_clusters = SymmetricClusterSet(lattice, square_symmetries)
                        clusters_from_clusters!(sym_clusters, trans_clusters)

                        @test length(sym_clusters) == 533
                end

                @testset "Isomorphic clustering (Square)" begin
                        lattice = SiteExpansionLattice(2, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)

                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)

                        order2 = filter(c -> length(c) == 2, collect(iso_clusters))
                        @test length(order2) == 1
                        @test order2[1].lc > 1
                end

        end # Clusters

        # -------------------------------------------------------------------
        @testset verbose = true "Expansions" begin

                @testset "Square Lattice pipeline" begin
                        m_order = 3
                        lattice = SiteExpansionLattice(m_order, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)

                        expansion = Expansion(iso_clusters, lattice, m_order)
                        summation!(expansion, m_order)

                        @test expansion.weights == [1.0 -4.0 6.0; 0.0 2.0 -12.0; 0.0 0.0 6.0]
                end

                @testset "Kagome Lattice pipeline" begin
                        m_order = 3
                        lattice = SiteExpansionLattice(m_order, kagome_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)

                        expansion = Expansion(iso_clusters, lattice, m_order)
                        summation!(expansion, m_order)

                        @test isapprox(expansion.weights, [1.0 -4 6.0; 0.0 2 -10.0; 0.0 0.0 0.6666666666666666; 0.0 0.0 4.0])
                end

                @testset "write_to_json" begin
                        m_order = 2
                        lattice = SiteExpansionLattice(m_order, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice, m_order)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, iso_clusters, path)

                                data = JSON3.read(read(path, String))

                                @test length(data) == length(iso_clusters)

                                for entry in data
                                        @test haskey(entry, :cluster_id)
                                        @test haskey(entry, :n_sites)
                                        @test haskey(entry, :coordinates)
                                        @test haskey(entry, :site_colors)
                                        @test haskey(entry, :bonds)
                                        @test haskey(entry, :weights)
                                end

                                order1 = filter(e -> e[:n_sites] == 1, collect(data))
                                @test length(order1) == 1
                                @test length(order1[1][:bonds]) == 0

                                order2 = filter(e -> e[:n_sites] == 2, collect(data))
                                @test length(order2) == 1
                                @test length(order2[1][:bonds]) == 1
                                @test order2[1][:bonds][1][end] == 1
                        end
                end

        end # Expansions

end # LINCEGE.jl
