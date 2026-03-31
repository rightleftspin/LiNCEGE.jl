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
