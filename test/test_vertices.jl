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
