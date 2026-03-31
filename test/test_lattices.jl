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
