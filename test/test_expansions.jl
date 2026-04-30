@testset verbose = true "Expansions" begin

        @testset "Square Lattice pipeline" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, square_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set(values(weights(expansion, 1))) == Set((1.0))
                @test Set(values(weights(expansion, 2))) == Set((-4.0, 2.0))
                @test Set(values(weights(expansion, 3))) == Set((6.0, -12.0, 6.0))

        end

        @testset "Kagome Lattice pipeline" begin
                m_order = 3
                lattice = SiteExpansionLattice(m_order, kagome_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((-4.0, 2.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((6.0, -10.0, round(2 / 3, digits=6), 4.0))
        end

        @testset "Pyrochlore Lattice Unit Cell pipeline" begin
                m_order = 3
                lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((0.25, -1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((-3.0, 1.5, 0.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 4))]) == Set((16.5, -27.0, 10.5, 1.0, 1.0, 0.0))
        end

        @testset "Square Lattice Cluster pipeline" begin
                m_order = 3
                lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)
                trans_clusters = TranslationClusterSet(lattice)
                clusters_from_lattice!(trans_clusters, lattice)
                iso_clusters = IsomorphicClusterSet(lattice)
                clusters_from_clusters!(iso_clusters, trans_clusters)

                expansion = Expansion(iso_clusters, lattice)
                summation!(expansion, m_order)

                @test Set([round(v, digits=6) for v in values(weights(expansion, 1))]) == Set((1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 2))]) == Set((0.5, -2.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 3))]) == Set((-2.0, 1.0, 1.0))
                @test Set([round(v, digits=6) for v in values(weights(expansion, 4))]) == Set((3.0, -6.0, 2.0, 1.0, 0.0))
        end

        @testset "write_to_json" begin
                expected_keys = ["bonds", "cluster_hash", "coordinates", "lattice_constant",
                        "n_sites", "order", "site_colors", "subgraphs", "weights"]

                @testset "SiteExpansionLattice" begin
                        m_order = 2
                        lattice = SiteExpansionLattice(m_order, square_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["n_sites"] == 1, data)
                                @test length(order1) == 1
                                @test length(order1[1]["bonds"]) == 0

                                order2 = filter(e -> e["n_sites"] == 2, data)
                                @test length(order2) == 1
                                @test length(order2[1]["bonds"]) == 1
                                @test order2[1]["bonds"][1][end] == 1
                        end
                end

                @testset "StrongClusterExpansionLattice" begin
                        m_order = 2
                        lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters) + n_unique_sites(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["order"] == 1, data)
                                @test !isempty(order1)
                                @test all(e -> e["n_sites"] == 1, order1)
                        end
                end

                @testset "WeakClusterExpansionLattice" begin
                        m_order = 2
                        lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)
                        trans_clusters = TranslationClusterSet(lattice)
                        clusters_from_lattice!(trans_clusters, lattice)
                        iso_clusters = IsomorphicClusterSet(lattice)
                        clusters_from_clusters!(iso_clusters, trans_clusters)
                        expansion = Expansion(iso_clusters, lattice)
                        summation!(expansion, m_order)

                        mktempdir() do dir
                                path = joinpath(dir, "expansion.json")
                                write_to_json(expansion, lattice, path)
                                data = JSON.parse(read(path, String))

                                @test length(data) == length(iso_clusters) + n_unique_sites(iso_clusters)

                                for entry in data
                                        @test sort(collect(keys(entry))) == expected_keys
                                end

                                order1 = filter(e -> e["order"] == 1, data)
                                @test !isempty(order1)
                                @test all(e -> e["n_sites"] == 1, order1)
                        end
                end
        end

end # Expansions
