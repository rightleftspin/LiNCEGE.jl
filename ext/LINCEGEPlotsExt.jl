module LINCEGEPlotsExt

using Plots
using ColorSchemes
import LINCEGE: UnitCells

function UnitCells.image_unit_cell(unit_cell::UnitCells.UnitCell)
        dim = UnitCells.dimension(unit_cell)
        try
                if dim == 2
                        return _image_unit_cell_2d(unit_cell)
                elseif dim == 3
                        gr()
                        return _image_unit_cell_3d(unit_cell)
                else
                        error("Only 2D and 3D unit cells can be visualized. Got dimension: $dim")
                end
        catch e
                @warn "Visualization failed: $e"
                return 1
        end
end

function _image_unit_cell_2d(unit_cell::UnitCells.UnitCell)
        sitesInfo = []

        maxXValue = 0.1
        minXValue = -0.1
        maxYValue = 0.1
        minYValue = -0.1

        pv1 = unit_cell.primitive_vectors[:, 1]
        pv2 = unit_cell.primitive_vectors[:, 2]

        unitcellLines = [
                [[0, 0], pv1],
                [[0, 0], pv2],
                [pv2, pv1 .+ pv2],
                [pv1, pv1 .+ pv2]
        ]

        unitCellPlot = plot(legend=:best)
        firstLine = true

        for i in 1:UnitCells.basis_size(unit_cell)
                basis_pos = unit_cell.basis[:, i]
                push!(sitesInfo, [basis_pos[1], basis_pos[2], unit_cell.site_colors[i]])
        end

        for bond in unit_cell.bonds
                prX = bond.direction[1]
                prY = bond.direction[2]

                basis1 = unit_cell.basis[:, bond.site1]
                basis2 = unit_cell.basis[:, bond.site2]

                site1X = basis1[1]
                site1Y = basis1[2]
                site2X = basis2[1] + prX * pv1[1] + prY * pv2[1]
                site2Y = basis2[2] + prX * pv1[2] + prY * pv2[2]

                maxXValue = max(maxXValue, site1X, site2X)
                minXValue = min(minXValue, site1X, site2X)
                maxYValue = max(maxYValue, site1Y, site2Y)
                minYValue = min(minYValue, site1Y, site2Y)

                my_palette = colorschemes[:seaborn_bright]
                plot!(unitCellPlot, [site1X, site2X], [site1Y, site2Y],
                        label="", color=my_palette[bond.bond_type], lw=2, dpi=1000)
        end

        minXValue = min(minXValue, pv1[1], pv2[1])
        maxXValue = max(maxXValue, pv1[1] + pv2[1], pv1[1], pv2[1])
        minYValue = min(minYValue, pv1[2], pv2[2])
        maxYValue = max(maxYValue, pv1[2] + pv2[2])

        xPositions = getindex.(sitesInfo, 1)
        yPositions = getindex.(sitesInfo, 2)
        colors = getindex.(sitesInfo, 3)
        paletteAtoms = colorschemes[:Accent_8]

        for c in unique(colors)
                mask = colors .== c
                colorIndex = Int(c)
                plot!(unitCellPlot, xPositions[mask], yPositions[mask],
                        seriestype=:scatter,
                        aspect_ratio=1,
                        xlims=(minXValue * 1.1, maxXValue * 1.1),
                        ylims=(minYValue * 1.1, maxYValue * 1.1),
                        color=paletteAtoms[(colorIndex+1)],
                        label="Type $colorIndex",
                        markerstrokewidth=2,
                        markersize=6,
                        dpi=1000)
        end

        for unitCellLine in unitcellLines
                xPoints = [unitCellLine[1][1], unitCellLine[2][1]]
                yPoints = [unitCellLine[1][2], unitCellLine[2][2]]
                plot!(unitCellPlot, xPoints, yPoints, linestyle=:dash, color="black",
                        label=firstLine ? "Unit Cell" : "")
                firstLine = false
        end

        display(unitCellPlot)
        return unitCellPlot
end

function _image_unit_cell_3d(unit_cell::UnitCells.UnitCell)
        sitesInfo = []
        maxXValue = 0.1; minXValue = -0.1
        maxYValue = 0.1; minYValue = -0.1
        maxZValue = 0.1; minZValue = -0.1

        pv1 = unit_cell.primitive_vectors[:, 1]
        pv2 = unit_cell.primitive_vectors[:, 2]
        pv3 = unit_cell.primitive_vectors[:, 3]

        unitcellLines = [
                [[0,0,0], pv1], [[0,0,0], pv2], [pv1, pv1.+pv2], [pv2, pv1.+pv2],
                [pv3, pv3.+pv1], [pv3, pv3.+pv2], [pv3.+pv1, pv3.+pv1.+pv2],
                [pv3.+pv2, pv3.+pv1.+pv2], [[0,0,0], pv3], [pv1, pv1.+pv3],
                [pv2, pv2.+pv3], [pv1.+pv2, pv1.+pv2.+pv3]
        ]

        unitCellPlot = plot3d(legend=:outertopright)
        firstLine = true

        for i in 1:UnitCells.basis_size(unit_cell)
                basis_pos = unit_cell.basis[:, i]
                push!(sitesInfo, [basis_pos[1], basis_pos[2], basis_pos[3], unit_cell.site_colors[i]])
        end

        for bond in unit_cell.bonds
                basis1 = unit_cell.basis[:, bond.site1]
                basis2 = unit_cell.basis[:, bond.site2]
                prX, prY, prZ = bond.direction[1], bond.direction[2], bond.direction[3]

                site1X, site1Y, site1Z = basis1[1], basis1[2], basis1[3]
                site2X = basis2[1] + prX*pv1[1] + prY*pv2[1] + prZ*pv3[1]
                site2Y = basis2[2] + prX*pv1[2] + prY*pv2[2] + prZ*pv3[2]
                site2Z = basis2[3] + prX*pv1[3] + prY*pv2[3] + prZ*pv3[3]

                maxXValue = max(maxXValue, site1X, site2X); minXValue = min(minXValue, site1X, site2X)
                maxYValue = max(maxYValue, site1Y, site2Y); minYValue = min(minYValue, site1Y, site2Y)
                maxZValue = max(maxZValue, site1Z, site2Z); minZValue = min(minZValue, site1Z, site2Z)

                my_palette = colorschemes[:seaborn_bright]
                plot3d!(unitCellPlot, [site1X, site2X], [site1Y, site2Y], [site1Z, site2Z],
                        label="", color=my_palette[bond.bond_type], lw=2)
        end

        for i in 1:3
                pv = unit_cell.primitive_vectors[:, i]
                minXValue = min(minXValue, pv[1]); maxXValue = max(maxXValue, pv[1])
                minYValue = min(minYValue, pv[2]); maxYValue = max(maxYValue, pv[2])
                minZValue = min(minZValue, pv[3]); maxZValue = max(maxZValue, pv[3])
        end

        for l in unitcellLines
                plot3d!(unitCellPlot, [l[1][1], l[2][1]], [l[1][2], l[2][2]], [l[1][3], l[2][3]],
                        linestyle=:dash, color="black", label=firstLine ? "Unit Cell" : "")
                firstLine = false
        end

        plot3d!(unitCellPlot,
                xlims=(minXValue*1.1, maxXValue*1.1),
                ylims=(minYValue*1.1, maxYValue*1.1),
                zlims=(minZValue*1.1, maxZValue*1.1),
                aspect_ratio=:equal, xlabel="X", ylabel="Y", zlabel="Z")

        xPositions = getindex.(sitesInfo, 1)
        yPositions = getindex.(sitesInfo, 2)
        zPositions = getindex.(sitesInfo, 3)
        colors = getindex.(sitesInfo, 4)
        paletteAtoms = colorschemes[:Accent_8]

        for c in unique(colors)
                mask = colors .== c
                colorIndex = Int(c)
                plot3d!(unitCellPlot, xPositions[mask], yPositions[mask], zPositions[mask],
                        seriestype=:scatter3d, markershape=:circle, aspect_ratio=:equal,
                        color=paletteAtoms[colorIndex], label="Type $colorIndex",
                        markerstrokewidth=2, markersize=4)
        end

        display(unitCellPlot)
        return unitCellPlot
end

end # module
