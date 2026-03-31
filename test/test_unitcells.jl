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
