using Base: find_all_in_cache_path
using Test
using Plots

@testset verbose = true "LINCEGE.jl" begin

        @testset verbose = true "Vertices" begin
                import LINCEGE:
                        Vertices.LatticeVertices,
                        Vertices.ExpansionVertices

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
        end

        @testset verbose = true "UnitCells" begin
                import LINCEGE:
                        UnitCells.UnitCell,
                        UnitCells.Bond,
                        UnitCells.dimension,
                        UnitCells.basis_size,
                        UnitCells.shift_unit_cell,
                        UnitCells.find_possible_neighbors,
                        UnitCells.image_unit_cell

                @testset "Square Lattice" begin
                        basis = [[0.0, 0.0]]
                        primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
                        bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
                        colors = [1]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)

                        @test basis_size(unit_cell) == length(basis)
                        @test dimension(unit_cell) == length(primitive_vectors)
                        @test shift_unit_cell(unit_cell, [1, 0, 1]) == [1.0, 0.0]
                        @test shift_unit_cell(unit_cell, [2, -2, 1]) == [2.0, -2.0]
                        @test shift_unit_cell(unit_cell, [1 2 3; 0 -2 1; 1 1 1]) == [1.0 2.0 3.0; 0.0 -2.0 1.0]
                        @test find_possible_neighbors(unit_cell, [0, 0, 1]) == [[1, 0, 1], [0, 1, 1]]
                        @test image_unit_cell(unit_cell) isa Plots.Plot{Plots.GRBackend}
                end

                @testset "Kagome Lattice" begin
                        basis = [[0.0, 0.0], [1.0, 0.0], [1 / 2, sqrt(3) / 2]]
                        colors = [1, 1, 1]

                        primitive_vectors = [[2, 0], [1, sqrt(3)]]
                        bonds = [Bond(1, 2, [0, 0], 1), Bond(1, 3, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)

                        @test image_unit_cell(unit_cell) isa Plots.Plot{Plots.GRBackend}
                        #The function displayes but also returns the plot. In order to save need to use a separate line
                        #savefig(image_unit_cell(unit_cell),"kagomeUnitCell.png")
                end

                @testset "One Fifth Depleted" begin
                        primitive_vectors = [[2.0, 1.0], [-1.0, 2.0]]
                        basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0], [1.0, 1.0]]
                        bonds = [Bond(1, 2, [0, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 4, [0, 0], 1), Bond(3, 4, [0, 1], 1), Bond(2, 3, [0, 0], 2), Bond(1, 4, [-1, 0], 2)]
                        colors = [1, 1, 1, 1]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)
                        @test image_unit_cell(unit_cell) isa Plots.Plot{Plots.GRBackend}

                end
                @testset "Shastry Sutherland" begin
                        primitive_vectors = [[(sqrt(3) + 1) / sqrt(2), 0], [0, (sqrt(3) + 1) / sqrt(2)]]
                        basis = [[0.0, 0.0], [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))], [sqrt(2) / 2, sqrt(6) / 2], [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
                        bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
                                Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
                                Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)]
                        colors = [1, 2, 3, 4]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)
                        @test image_unit_cell(unit_cell) isa Plots.Plot{Plots.GRBackend}
                end
                @testset "3D Cube" begin
                        primitive_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                        basis = [[0.0, 0.0, 0.0]]
                        colors = [1]
                        bonds = [Bond(1, 1, [1, 0, 0], 1), Bond(1, 1, [0, 1, 0], 1), Bond(1, 1, [0, 0, 1], 1)]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)
                        @test image_unit_cell(unit_cell) isa Plots.Plot{Plots.GRBackend}
                end
        end

        @testset verbose = true "Lattices" begin
                import LINCEGE:
                        Lattices.AbstractLattice,
                        Lattices.SiteExpansionLattice,
                        Lattices.centers,
                        Lattices.max_order,
                        Lattices.neighbors,
                        Lattices.n_unique_sites

                @testset verbose = true "InfiniteLattice" begin
                        @testset verbose = true "SiteExpansionLattice" begin
                                @testset "Square Lattice" begin
                                        basis = [[0.0, 0.0]]
                                        primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
                                        bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
                                        colors = [1]
                                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)

                                        m_order = 2
                                        lattice = SiteExpansionLattice(m_order, unit_cell)

                                        @test centers(lattice) == ExpansionVertices([13])
                                        @test max_order(lattice) == UInt8(2)
                                        @test n_unique_sites(lattice) == 1
                                        @test neighbors(lattice, centers(lattice)) == ExpansionVertices([12, 14, 8, 18])
                                end
                                @testset "Kagome Lattice" begin
                                        basis = [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]]
                                        colors = [1, 1, 1]
                                        primitive_vectors = [[2, 0], [1, sqrt(3)]]
                                        bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1), Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [0, 1], 1)]
                                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)

                                        m_order = 2
                                        lattice = SiteExpansionLattice(m_order, unit_cell)

                                        @test centers(lattice) == ExpansionVertices([13, 38, 63])
                                        @test max_order(lattice) == UInt8(2)
                                        @test n_unique_sites(lattice) == 3
                                        @test neighbors(lattice, ExpansionVertices([38])) == ExpansionVertices([13, 63, 14, 68])
                                end
                        end
                end
        end
        @testset verbose = true "Clusters" begin
                import LINCEGE:
                        Clusters.ClusterSet,
                        Clusters.IsomorphicClusterSet,
                        Clusters.clusters_from_lattice!,
                        Clusters.clusters_from_clusters!,
                        Expansions.SiteExpansion,
                        Expansions.summation!

                @testset "Translation Clusters Square" begin
                        basis = [[0.0, 0.0]]
                        primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
                        bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
                        colors = [1]
                        unit_cell = UnitCell(basis, primitive_vectors, bonds, colors)

                        m_order = 3
                        lattice = SiteExpansionLattice(m_order, unit_cell)
                        translation_clusters = ClusterSet(lattice)
                        clusters_from_lattice!(translation_clusters, lattice)

                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, translation_clusters)

                        expansion = SiteExpansion(iso_clusters, lattice, m_order)
                        summation!(expansion, m_order)
                        @test expansion.weights == [1.0 -4.0 6.0; 0.0 2.0 -12.0; 0.0 0.0 6.0]
                end

        end
end
