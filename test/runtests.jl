using Test

@testset verbose = true "LINCEGE.jl" begin
    import LINCEGE:
        Lattices.SiteExpansionLattice,
        Lattices.StrongClusterExpansionLattice,
        Lattices.n_unique_sites,
        ClusterCollections.TranslationClusters,
        ClusterCollections.IsomorphicClusters,
        ClusterExpansions.SiteExpansion,
        ClusterExpansions.summation!

    basis = [[0.0, 0.0]]
    primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
    colors = [1]
    neighbor_distances = [1.0]
    max_order_square = 3

    square_lattice = SiteExpansionLattice(
        max_order_square,
        basis,
        primitive_vectors,
        colors,
        neighbor_distances,
    )
    translation_clusters_square = TranslationClusters(square_lattice)
    isomorphic_clusters_square = IsomorphicClusters(translation_clusters_square, square_lattice)
    square_expansion = SiteExpansion(isomorphic_clusters_square, square_lattice)
    summation!(square_expansion, isomorphic_clusters_square, square_lattice)

    basis = [[0.0, 0.0]]
    primitive_vectors = [[0.5, sqrt(3) / 2], [1, 0]]
    colors = [1]
    neighbor_distances = [1.0]
    max_order_triangular = 3

    triangular_lattice = SiteExpansionLattice(
        max_order_triangular,
        basis,
        primitive_vectors,
        colors,
        neighbor_distances,
    )
    translation_clusters_triangular = TranslationClusters(triangular_lattice)
    isomorphic_clusters_triangular = IsomorphicClusters(translation_clusters_triangular, triangular_lattice)
    triangular_expansion = SiteExpansion(isomorphic_clusters_triangular, triangular_lattice)
    summation!(triangular_expansion, isomorphic_clusters_triangular, triangular_lattice)

    basis = [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]]
    primitive_vectors = [[2, 0], [1, sqrt(3)]]
    colors = [1, 1, 1]
    neighbor_distances = [1.0]
    max_order_kagome = 3

    kagome_lattice = SiteExpansionLattice(
        max_order_kagome,
        basis,
        primitive_vectors,
        colors,
        neighbor_distances,
    )
    translation_clusters_kagome = TranslationClusters(kagome_lattice)
    isomorphic_clusters_kagome = IsomorphicClusters(translation_clusters_kagome, kagome_lattice)
    kagome_expansion = SiteExpansion(isomorphic_clusters_kagome, kagome_lattice)
    summation!(kagome_expansion, isomorphic_clusters_kagome, kagome_lattice)

    tiling_units = [[
        [1 / 2, 1 / 2, 1 / 2],
        [1 / 2, -1 / 2, -1 / 2],
        [-1 / 2, 1 / 2, -1 / 2],
        [-1 / 2, -1 / 2, 1 / 2],
    ]]
    primitive_vectors = [[2, 2, 0], [2, 0, 2], [0, 2, 2]]
    colors = [[1, 1, 1, 1]]
    neighbor_distances = [sqrt(2)]
    max_order_pyrochlore_unit = 7

    pyrochlore_unit_lattice = StrongClusterExpansionLattice(
        max_order_pyrochlore_unit,
        tiling_units,
        primitive_vectors,
        colors,
        neighbor_distances,
    )
    translation_clusters_pyrochlore_unit = TranslationClusters(pyrochlore_unit_lattice)
    isomorphic_clusters_pyrochlore_unit = IsomorphicClusters(translation_clusters_pyrochlore_unit, pyrochlore_unit_lattice)
    #pyrochlore_unit_expansion = ClusterExpansion(isomorphic_clusters_pyrochlore_unit, pyrochlore_unit_lattice)
    #summation!(pyrochlore_unit_expansion, isomorphic_clusters_pyrochlore_unit, pyrochlore_unit_lattice)

    @testset "Vertices" begin
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

    @testset "Lattices" begin
        import LINCEGE:
            Vertices.LatticeVertices,
            Lattices.SiteExpansionLattice,
            Lattices.centers,
            Lattices.max_order,
            Lattices.neighbors

        @testset "SiteExpansionLattice" begin
            @testset "Square Lattice" begin
                size_square = max_order_square * 2 + 1
                center = size_square * fld(size_square, 2) + round(Int, size_square / 2, RoundUp)
                @test centers(square_lattice) == LatticeVertices(center)
                @test length(neighbors(square_lattice, LatticeVertices(center))) == 4
                @test n_unique_sites(square_lattice) == 1
            end

            @testset "Triangular Lattice" begin
                size_triangular = max_order_triangular * 2 + 1
                center = size_triangular * fld(size_triangular, 2) + round(Int, size_triangular / 2, RoundUp)
                @test centers(triangular_lattice) == LatticeVertices(center)
                @test length(neighbors(triangular_lattice, LatticeVertices(center))) == 6
                @test n_unique_sites(triangular_lattice) == 1
            end

            @testset "Kagome Lattice" begin
                size_kagome = max_order_kagome * 2 + 1
                center = size_kagome * fld(size_kagome, 2) + round(Int, size_kagome / 2, RoundUp)
                kagome_centers = [0, size_kagome^2, 2 * size_kagome^2] .+ center
                @test centers(kagome_lattice) == LatticeVertices(kagome_centers)
                @test all([length(neighbors(kagome_lattice, LatticeVertices(v))) == 4 for v in centers(kagome_lattice)])
                @test n_unique_sites(kagome_lattice) == 3
            end
        end
    end

    @testset "StrongClusterExpansionlattice" begin
        import LINCEGE:
            Vertices.LatticeVertices,
            Vertices.ExpansionVertices,
            Lattices.centers,
            Lattices.max_order,
            Lattices.neighbors

        @testset "Pyrochlore Lattice Unit Cell Expansion" begin
            size_pyrochlore_unit = max_order_pyrochlore_unit * 2 + 1
            center = round(Int, (size_pyrochlore_unit^3) / 2, RoundUp)
            @test centers(pyrochlore_unit_lattice) == ExpansionVertices(center)
            @test all([length(neighbors(pyrochlore_unit_lattice, ExpansionVertices(v))) == 12 for v in centers(pyrochlore_unit_lattice)])
            @test n_unique_sites(pyrochlore_unit_lattice) == 4
        end
    end

    @testset "ClusterCollections" begin
        import LINCEGE:
            Clusters.lattice_constant
        @testset "TranslationClusters" begin
            @testset "Square Lattice" begin
                # Counts according to Rigol, Bryant and Singh 2007
                sum_lc_clusters = [1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446]
                @test length(translation_clusters_square) == sum(sum_lc_clusters[1:max_order_square])
            end
            @testset "Triangular Lattice" begin
                # Counts according to Rigol, Bryant and Singh 2007
                sum_lc_clusters = [1, 3, 11, 44, 186, 814, 3652, 16689, 77359, 362671]
                @test length(translation_clusters_triangular) == sum(sum_lc_clusters[1:max_order_triangular])
            end
            @testset "kagome lattice" begin
                # counts according to rigol, bryant and singh 2007
                sum_lc_clusters = [1, 2, 14 // 3, 12, 33, 281 // 3, 272, 805, 2420, 7358]
                @test length(translation_clusters_kagome) == 3 * sum(sum_lc_clusters[1:max_order_kagome])
            end
            @testset "pyrochlore unit lattice" begin
                # counts according to Schefer et al. 2020
                sum_lc_clusters = [4, 1, 6, 50, 475, 4881, 52835, 593382, 6849415] .// 4
                @test length(translation_clusters_pyrochlore_unit) == 4 * sum(sum_lc_clusters[1:max_order_pyrochlore_unit+1])
            end
        end
        @testset "IsomorphicClusters" begin
            @testset "Square Lattice" begin
                # Counts according to Rigol, Bryant and Singh 2007
                num_topo_clusters = [1, 1, 1, 3, 4, 10, 19, 51, 112, 300]
                sum_lc_clusters = [1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446]
                @test length(isomorphic_clusters_square) == sum(num_topo_clusters[1:max_order_square])
                @test sum([lattice_constant(cluster) for (ghash, cluster) in isomorphic_clusters_square]) == sum(sum_lc_clusters[1:max_order_square])
            end
            @testset "Triangular Lattice" begin
                # Counts according to Rigol, Bryant and Singh 2007
                num_topo_clusters = [1, 1, 2, 4, 8, 22, 54, 156, 457, 1424]
                sum_lc_clusters = [1, 3, 11, 44, 186, 814, 3652, 16689, 77359, 362671]
                @test length(isomorphic_clusters_triangular) == sum(num_topo_clusters[1:max_order_triangular])
                @test sum([lattice_constant(cluster) for (ghash, cluster) in isomorphic_clusters_triangular]) == sum(sum_lc_clusters[1:max_order_triangular])
            end
            @testset "Kagome Lattice" begin
                # Counts according to Rigol, Bryant and Singh 2007
                num_topo_clusters = [1, 1, 2, 2, 4, 7, 12, 22, 45, 88]
                sum_lc_clusters = [1, 2, 14 // 3, 12, 33, 281 // 3, 272, 805, 2420, 7358]
                @test length(isomorphic_clusters_kagome) == sum(num_topo_clusters[1:max_order_kagome])
                @test sum([lattice_constant(cluster) for (ghash, cluster) in isomorphic_clusters_kagome]) == 3 * sum(sum_lc_clusters[1:max_order_kagome])
            end
            @testset "Pyrochlore Unit Lattice" begin
                # counts according to Schefer et al. 2020
                num_topo_clusters = [1, 1, 1, 3, 8, 25, 100, 466, 2473]
                sum_lc_clusters = [4, 1, 6, 50, 475, 4881, 52835, 593382, 6849415] .// 4
                @test length(isomorphic_clusters_pyrochlore_unit) == sum(num_topo_clusters[1:max_order_pyrochlore_unit+1])
                @test sum([lattice_constant(cluster) for (ghash, cluster) in isomorphic_clusters_pyrochlore_unit]) == 4 * sum(sum_lc_clusters[1:max_order_pyrochlore_unit+1])
            end
        end
    end
end
