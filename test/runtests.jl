using NLCE, Test

@testset verbose = true "NLCE" begin

    @testset "Square Lattice Site Expansion" begin

        square_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[1, 0], [0, 1]])
        # Counts according to Rigol, Bryant and Singh 2007

        num_topo_clusters_square = [1, 1, 1, 3, 4, 10, 19, 51, 112, 300]
        sum_lc_clusters_square = [1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446]

        neighbors = [1]
        max_order = 6

        square_nlce_bundle = NLCE.SiteExpansionBundle(
            square_lattice["Basis"],
            square_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            square_nlce_bundle,
            length(NLCE.start(square_nlce_bundle)),
        )

        square_lattice_cluster_info = NLCE.lattice_constants!(
            square_nlce_bundle,
            length(NLCE.start(square_nlce_bundle)),
            t_i_clusters,
            super_verts,
        )

        @test sum([mult for (hash, (_, mult, _, _, _)) in square_lattice_cluster_info]) == sum(sum_lc_clusters_square[1:max_order])
        @test length(square_lattice_cluster_info) ==
              sum(num_topo_clusters_square[1:max_order])

        # TODO: Write tests for subcluster multiplicities

    end


    @testset "Triangular Lattice Site Expansion" begin

        triangular_lattice =
            Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[0.5, sqrt(3) / 2], [1, 0]])
        # Counts according to Rigol, Bryant and Singh 2007

        num_topo_clusters_triangular = [1, 1, 2, 4, 8, 22, 54, 156, 457, 1424]
        sum_lc_clusters_triangular = [1, 3, 11, 44, 186, 814, 3652, 16689, 77359, 362671]

        neighbors = [1]
        max_order = 6

        triangular_nlce_bundle = NLCE.SiteExpansionBundle(
            triangular_lattice["Basis"],
            triangular_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            triangular_nlce_bundle,
            length(NLCE.start(triangular_nlce_bundle)),
        )

        triangular_lattice_cluster_info = NLCE.lattice_constants!(
            triangular_nlce_bundle,
            length(NLCE.start(triangular_nlce_bundle)),
            t_i_clusters,
            super_verts,
        )

        @test sum([
            mult for (hash, (_, mult, _, _, _)) in triangular_lattice_cluster_info
        ]) == sum(sum_lc_clusters_triangular[1:max_order])
        @test length(triangular_lattice_cluster_info) ==
              sum(num_topo_clusters_triangular[1:max_order])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Kagome Lattice Site Expansion" begin

        kagome_lattice = Dict(
            "Basis" => [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]],
            "Primitive Vectors" => [[2, 0], [1, sqrt(3)]],
        )
        # Counts according to Rigol, Bryant and Singh 2007

        num_topo_clusters_kagome = [1, 1, 2, 2, 4, 7, 12, 22, 45, 88]
        sum_lc_clusters_kagome = [1, 2, 14 // 3, 12, 33, 281 // 3, 272, 805, 2420, 7358]

        neighbors = [1]
        max_order = 6

        kagome_nlce_bundle = NLCE.SiteExpansionBundle(
            kagome_lattice["Basis"],
            kagome_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            kagome_nlce_bundle,
            length(NLCE.start(kagome_nlce_bundle)),
        )

        kagome_lattice_cluster_info = NLCE.lattice_constants!(
            kagome_nlce_bundle,
            length(NLCE.start(kagome_nlce_bundle)),
            t_i_clusters,
            super_verts,
        )


        @test sum([mult for (hash, (_, mult, _, _, _)) in kagome_lattice_cluster_info]) == sum(sum_lc_clusters_kagome[1:max_order])
        @test length(kagome_lattice_cluster_info) ==
              sum(num_topo_clusters_kagome[1:max_order])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Square Lattice Cluster Expansion" begin

        square_lattice = Dict(
            "Expansion Basis" => [[1 / 2, 1 / 2]],
            "Struct Per Basis" =>
                [[[-1 / 2, -1 / 2], [-1 / 2, 1 / 2], [1 / 2, -1 / 2], [1 / 2, 1 / 2]]],
            "Expansion Labels" => [[1, 2, 2, 1]],
            "Expansion Primitive Vectors" => [[1, 1], [1, -1]],
            "Expansion Neighbors" => [sqrt(2)],
        )

        # Counts according to Rigol, Bryant and Singh 2007
        num_topo_clusters_square = [1, 1, 1, 2, 5, 11]
        sum_lc_clusters_square = [1, 1 // 2, 1, 3, 19 // 2, 63 // 2]

        # Order starts a 0 for the single site, then goes up from there
        neighbors = [1]
        max_order = 5

        square_nlce_bundle = NLCE.WeakClusterExpansionBundle(
            square_lattice["Expansion Basis"],
            square_lattice["Struct Per Basis"],
            square_lattice["Expansion Labels"],
            square_lattice["Expansion Primitive Vectors"],
            square_lattice["Expansion Neighbors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            square_nlce_bundle,
            (length(unique(Iterators.flatten(square_lattice["Expansion Labels"])))),
            single_site = true,
        )

        square_lattice_cluster_info = NLCE.lattice_constants!(
            square_nlce_bundle,
            (length(unique(Iterators.flatten(square_lattice["Expansion Labels"])))),
            t_i_clusters,
            super_verts,
        )

        # Because counting starts at 0, we have to check from the first count to the max_order + 1 count
        @test sum([mult for (hash, (_, mult, _, _, _)) in square_lattice_cluster_info]) == sum(sum_lc_clusters_square[1:(max_order+1)])
        @test length(square_lattice_cluster_info) ==
              sum(num_topo_clusters_square[1:(max_order+1)])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Triangular Lattice Cluster Expansion" begin

        triangular_lattice = Dict(
            "Expansion Basis" => [[1 / 2, sqrt(3) / 4]],
            "Struct Per Basis" =>
                [[[-1 / 2, -sqrt(3) / 4], [1 / 2, -sqrt(3) / 4], [0, sqrt(3) / 4]]],
            "Expansion Labels" => [[1, 1, 1]],
            "Expansion Primitive Vectors" => [[1, 0], [1 / 2, sqrt(3) / 2]],
            "Expansion Neighbors" => [1],
        )

        # Counts according to Rigol, Bryant and Singh 2007
        num_topo_clusters_triangular = [1, 1, 1, 3, 5, 12, 35, 98, 299]
        sum_lc_clusters_triangular =
            [1, 1 // 3, 1, 11 // 3, 44 // 3, 62, 814 // 3, 3652 // 3, 5563]

        # Order starts a 0 for the single site, then goes up from there
        neighbors = [1]
        max_order = 5

        triangular_nlce_bundle = NLCE.WeakClusterExpansionBundle(
            triangular_lattice["Expansion Basis"],
            triangular_lattice["Struct Per Basis"],
            triangular_lattice["Expansion Labels"],
            triangular_lattice["Expansion Primitive Vectors"],
            triangular_lattice["Expansion Neighbors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts =
            NLCE.initial_clusters(triangular_nlce_bundle, 3, single_site = true)

        triangular_lattice_cluster_info =
            NLCE.lattice_constants!(triangular_nlce_bundle, 3, t_i_clusters, super_verts)


        # Because counting starts at 0, we have to check from the first count to the max_order + 1 count
        @test sum([
            mult for (hash, (_, mult, _, _, _)) in triangular_lattice_cluster_info
        ]) == sum(sum_lc_clusters_triangular[1:(max_order+1)])
        @test length(triangular_lattice_cluster_info) ==
              sum(num_topo_clusters_triangular[1:(max_order+1)])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Kagome Lattice Cluster Expansion" begin

        kagome_lattice = Dict(
            "Expansion Basis" => [[0, sqrt(3) / 3], [0, -sqrt(3) / 3]],
            "Struct Per Basis" => [
                [[0, -sqrt(3) / 3], [1 / 2, sqrt(3) / 6], [-1 / 2, sqrt(3) / 6]],
                [[0, sqrt(3) / 3], [1 / 2, -sqrt(3) / 6], [-1 / 2, -sqrt(3) / 6]],
            ],
            "Expansion Labels" => [[1, 2, 3], [1, 3, 2]],
            "Expansion Primitive Vectors" => [[1, sqrt(3)], [1, -sqrt(3)]],
            "Expansion Neighbors" => [2 * sqrt(3) / 3],
        )

        # Counts according to Rigol, Bryant and Singh 2007
        num_topo_clusters_kagome = [1, 1, 1, 1, 2, 2, 5, 7, 15]
        sum_lc_clusters_kagome = [1, 2 // 3, 1, 2, 14 // 3, 12, 94 // 3, 250 // 3, 225]

        # Order starts a 0 for the single site, then goes up from there
        neighbors = [1]
        max_order = 8

        kagome_nlce_bundle = NLCE.WeakClusterExpansionBundle(
            kagome_lattice["Expansion Basis"],
            kagome_lattice["Struct Per Basis"],
            kagome_lattice["Expansion Labels"],
            kagome_lattice["Expansion Primitive Vectors"],
            kagome_lattice["Expansion Neighbors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            kagome_nlce_bundle,
            (length(unique(Iterators.flatten(kagome_lattice["Expansion Labels"])))),
            single_site = true,
        )

        kagome_lattice_cluster_info = NLCE.lattice_constants!(
            kagome_nlce_bundle,
            (length(unique(Iterators.flatten(kagome_lattice["Expansion Labels"])))),
            t_i_clusters,
            super_verts,
        )

        # Because counting starts at 0, we have to check from the first count to the max_order + 1 count
        @test sum([mult for (hash, (_, mult, _, _, _)) in kagome_lattice_cluster_info]) == sum(sum_lc_clusters_kagome[1:(max_order+1)])
        @test length(kagome_lattice_cluster_info) ==
              sum(num_topo_clusters_kagome[1:(max_order+1)])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Pyrochlore Lattice Tetrahedra Expansion" begin

        pyrochlore_lattice = Dict(
            "Expansion Basis" => [[-1 / 2, -1 / 2, -1 / 2], [1 / 2, 1 / 2, 1 / 2]],
            "Struct Per Basis" => [
                [
                    [1 / 2, 1 / 2, 1 / 2],
                    [1 / 2, -1 / 2, -1 / 2],
                    [-1 / 2, 1 / 2, -1 / 2],
                    [-1 / 2, -1 / 2, 1 / 2],
                ],
                [
                    [-1 / 2, -1 / 2, -1 / 2],
                    [-1 / 2, 1 / 2, 1 / 2],
                    [1 / 2, -1 / 2, 1 / 2],
                    [1 / 2, 1 / 2, -1 / 2],
                ],
            ],
            "Expansion Labels" => [[1, 2, 3, 4], [1, 2, 3, 4]],
            "Expansion Primitive Vectors" => [[2, 2, 0], [2, 0, 2], [0, 2, 2]],
            "Expansion Neighbors" => [sqrt(3)],
        )

        # Counts according to Schäfer et al. 2020
        num_topo_clusters_pyrochlore = [1, 1, 1, 1, 2, 3, 6, 10, 24, 49]
        # Counts according to Applegate et al. 2012
        sum_lc_clusters_pyrochlore = [1, 1 // 2, 1, 3, 11]

        # Order starts a 0 for the single site, then goes up from there
        neighbors = [sqrt(2)]
        max_order = 4

        pyrochlore_nlce_bundle = NLCE.WeakClusterExpansionBundle(
            pyrochlore_lattice["Expansion Basis"],
            pyrochlore_lattice["Struct Per Basis"],
            pyrochlore_lattice["Expansion Labels"],
            pyrochlore_lattice["Expansion Primitive Vectors"],
            pyrochlore_lattice["Expansion Neighbors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            pyrochlore_nlce_bundle,
            (length(unique(Iterators.flatten(pyrochlore_lattice["Expansion Labels"])))),
            single_site = true,
        )

        pyrochlore_lattice_cluster_info = NLCE.lattice_constants!(
            pyrochlore_nlce_bundle,
            (length(unique(Iterators.flatten(pyrochlore_lattice["Expansion Labels"])))),
            t_i_clusters,
            super_verts,
        )

        # Because counting starts at 0, we have to check from the first count to the max_order + 1 count
        @test sum([
            mult for (hash, (_, mult, _, _, _)) in pyrochlore_lattice_cluster_info
        ]) == sum(sum_lc_clusters_pyrochlore[1:(max_order+1)])
        @test length(pyrochlore_lattice_cluster_info) ==
              sum(num_topo_clusters_pyrochlore[1:(max_order+1)])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Pyrochlore Lattice Unit Cell Expansion" begin

        pyrochlore_lattice = Dict(
            "Expansion Basis" => [[0, 0, 0]],
            "Struct Per Basis" => [[
                [1 / 2, 1 / 2, 1 / 2],
                [1 / 2, -1 / 2, -1 / 2],
                [-1 / 2, 1 / 2, -1 / 2],
                [-1 / 2, -1 / 2, 1 / 2],
            ]],
            "Expansion Labels" => [[1, 2, 3, 4]],
            "Expansion Primitive Vectors" => [[2, 2, 0], [2, 0, 2], [0, 2, 2]],
            "Expansion Neighbors" => [2 * sqrt(2)],
        )

        # Counts according to Schäfer et al. 2020
        num_topo_clusters_pyrochlore = [1, 1, 1, 3, 8, 25, 100, 466, 2473]

        # Order starts a 0 for the single site, then goes up from there
        neighbors = [sqrt(2)]
        max_order = 4

        pyrochlore_nlce_bundle = NLCE.StrongClusterExpansionBundle(
            pyrochlore_lattice["Expansion Basis"],
            pyrochlore_lattice["Struct Per Basis"],
            pyrochlore_lattice["Expansion Labels"],
            pyrochlore_lattice["Expansion Primitive Vectors"],
            pyrochlore_lattice["Expansion Neighbors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning,
        )

        t_i_clusters, super_verts = NLCE.initial_clusters(
            pyrochlore_nlce_bundle,
            (length(unique(Iterators.flatten(pyrochlore_lattice["Expansion Labels"])))),
            single_site = true,
        )

        pyrochlore_lattice_cluster_info = NLCE.lattice_constants!(
            pyrochlore_nlce_bundle,
            (length(unique(Iterators.flatten(pyrochlore_lattice["Expansion Labels"])))),
            t_i_clusters,
            super_verts,
        )

        # Because counting starts at 0, we have to check from the first count to the max_order + 1 count
        @test length(pyrochlore_lattice_cluster_info) ==
              sum(num_topo_clusters_pyrochlore[1:(max_order+1)])

        # TODO: Write tests for subcluster multiplicities

    end

end

@testset "Utility Functions" begin
    # Testing various utility functions
    test_adj_list = [[2, 3], [1, 4], [1, 4], [2, 3]]
    test_connections = [[1, 2], [2, 3], [4, 5], [5, 6]]

    test_adj_matrices = zeros(Int, 3, 6, 6)
    test_adj_matrices[1, :, :] = [
        1 1 0 0 0 0;
        1 1 1 0 0 0;
        0 1 1 0 0 0;
        0 0 0 1 1 0;
        0 0 0 1 1 1;
        0 0 0 0 1 1
    ]

    test_adj_matrices[2, :, :] = [
        1 1 0 0 0 0;
        3 1 1 0 0 0;
        0 3 1 0 0 0;
        0 0 0 1 1 0;
        0 0 0 3 1 1;
        0 0 0 0 3 1
    ]

    test_adj_matrices[3, :, :] = [
        0 1 0 0 0 0;
        1 0 2 0 0 0;
        0 2 0 0 0 0;
        0 0 0 0 3 0;
        0 0 0 3 0 4;
        0 0 0 0 4 0
    ]

    ans1 = zeros(Int, 3, 2, 2)
    ans1[1, :, :] = [
        1 1;
        1 1
    ]
    ans1[2, :, :] = [
        1 1;
        3 1
    ]
    ans1[3, :, :] = [
        0 1;
        1 0
    ]

    @test NLCE.reindex_adj_list(test_adj_list, [4]) == [[]]
    @test NLCE.reindex_adj_list(test_adj_list, [1, 2, 3]) == [[2, 3], [1], [1]]
    @test NLCE.reindex_adj_list(test_adj_list, [1, 2, 4]) == [[2], [1, 3], [2]]

    @test NLCE.reindex_connections(test_connections, [1], [1, 2]) == [[1, 2]]
    @test NLCE.reindex_connections(test_connections, [1, 2, 4], [1, 2, 3, 5, 6]) ==
          [[1, 2], [2, 3], [4, 5]]

    @test NLCE.reindex_adjacency_matrices(test_adj_matrices, [1], [1, 2]) == ans1

end
