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
