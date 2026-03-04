using Plots
using ColorSchemes
struct Bond
    site1::Int
    site2::Int
    direction::Vector{Int}
    bond_type::Int
end

neighbor_site(bond::Bond, coordinate::AbstractVector{Int}) = [coordinate[1:end-1] + bond.direction; bond.site2]

# coordinates are written [x1 x2 x3 ...; y1 y2 y3 ...; z1 z2 z3 ...]
struct UnitCell
    basis::Matrix{Float64}
    primitive_vectors::Matrix{Float64}
    site_colors::Vector{Int}
    bonds::Vector{Bond}
end

function UnitCell(basis::AbstractVector{<:AbstractVector{Float64}}, primitive_vectors::AbstractVector{<:AbstractVector{Float64}}, bonds::AbstractVector{Bond}, site_colors::AbstractVector{Int})

    @assert length(site_colors) == length(basis) "Number of site colors must match number of sites in the basis."
    @assert all(length.(primitive_vectors) .== length(primitive_vectors)) "Primitive vectors must form a square matrix"
    @assert all(length.(basis) .== length(basis[1])) "Basis vectors must all have the same dimension"


    UnitCell(hcat(basis...), hcat(primitive_vectors...), site_colors, bonds)
end

basis_size(unit_cell::UnitCell) = size(unit_cell.basis, 2)
dimension(unit_cell::UnitCell) = size(unit_cell.primitive_vectors, 2)
shift_unit_cell(unit_cell::UnitCell, shift_vector::AbstractVector{Int}) = unit_cell.primitive_vectors * shift_vector[1:end-1] + unit_cell.basis[:, shift_vector[end]]
shift_unit_cell(unit_cell::UnitCell, shift_vectors::AbstractMatrix{Int}) = stack(v -> shift_unit_cell(unit_cell, v), eachcol(shift_vectors))
find_possible_neighbors(unit_cell::UnitCell, coordinate::AbstractVector{Int}) = [neighbor_site(bond, coordinate) for bond in unit_cell.bonds if bond.site1 == coordinate[end]]
function image_unit_cell(unit_cell::UnitCell)
    dim = dimension(unit_cell)
    try
        if dim == 2
            image = image_unit_cell_2d(unit_cell)
            println(typeof(image))
            return image
        elseif dim == 3
            gr()  # Enable interactive 3D
            image = image_unit_cell_3d(unit_cell)
            println(typeof(image))
            return image
        else
            error("Only 2D and 3D unit cells can be visualized. Got dimension: $dim")
            return 1
        end
    catch e
        @warn "Visualization failed: $e"
        return 1
    end
end

function image_unit_cell_2d(unit_cell::UnitCell)
    sitesInfo = []

    maxXValue = 0.1
    minXValue = -0.1
    maxYValue = 0.1
    minYValue = -0.1

    # Get primitive vectors
    pv1 = unit_cell.primitive_vectors[:, 1]
    pv2 = unit_cell.primitive_vectors[:, 2]

    # Define unit cell boundary lines
    unitcellLines = [
        [[0, 0], pv1],
        [[0, 0], pv2],
        [pv2, pv1 .+ pv2],
        [pv1, pv1 .+ pv2]
    ]

    unitCellPlot = plot(legend=:best)
    firstLine = true


    # Plot basis sites
    for i in 1:basis_size(unit_cell)
        basis_pos = unit_cell.basis[:, i]
        locX = basis_pos[1]
        locY = basis_pos[2]
        color = unit_cell.site_colors[i]
        push!(sitesInfo, [locX, locY, color])
    end

    # Plot bonds
    for bond in unit_cell.bonds
        prX = bond.direction[1]
        prY = bond.direction[2]

        basis1 = unit_cell.basis[:, bond.site1]
        basis2 = unit_cell.basis[:, bond.site2]

        site1X = basis1[1]
        site1Y = basis1[2]
        site2X = basis2[1] + prX * pv1[1] + prY * pv2[1]
        site2Y = basis2[2] + prX * pv1[2] + prY * pv2[2]

        # Update min/max
        maxXValue = max(maxXValue, site1X, site2X)
        minXValue = min(minXValue, site1X, site2X)
        maxYValue = max(maxYValue, site1Y, site2Y)
        minYValue = min(minYValue, site1Y, site2Y)

        xData = [site1X, site2X]
        yData = [site1Y, site2Y]
        bondColor = bond.bond_type
        my_palette = colorschemes[:seaborn_bright]
        plot!(unitCellPlot, xData, yData, label="", color=my_palette[bondColor], lw=2, dpi=1000)
    end


    # Update bounds based on primitive vectors
    minXValue = min(minXValue, pv1[1], pv2[1])
    maxXValue = max(maxXValue, pv1[1] + pv2[1], pv1[1], pv2[1])
    minYValue = min(minYValue, pv1[2], pv2[2])
    maxYValue = max(maxYValue, pv1[2] + pv2[2])

    # Plot atoms
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

    # Plot unit cell edges
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

function image_unit_cell_3d(unit_cell::UnitCell)
    sitesInfo = []
    maxXValue = 0.1
    minXValue = -0.1
    maxYValue = 0.1
    minYValue = -0.1
    maxZValue = 0.1
    minZValue = -0.1

    # Get primitive vectors
    pv1 = unit_cell.primitive_vectors[:, 1]
    pv2 = unit_cell.primitive_vectors[:, 2]
    pv3 = unit_cell.primitive_vectors[:, 3]

    # Define the lines that form the 8 corners of the unit cell
    unitcellLines = [
        [[0, 0, 0], pv1],
        [[0, 0, 0], pv2],
        [pv1, pv1 .+ pv2],
        [pv2, pv1 .+ pv2],
        [pv3, pv3 .+ pv1],
        [pv3, pv3 .+ pv2],
        [pv3 .+ pv1, pv3 .+ pv1 .+ pv2],
        [pv3 .+ pv2, pv3 .+ pv1 .+ pv2],
        [[0, 0, 0], pv3],
        [pv1, pv1 .+ pv3],
        [pv2, pv2 .+ pv3],
        [pv1 .+ pv2, pv1 .+ pv2 .+ pv3]
    ]

    unitCellPlot = plot3d(legend=:outertopright)
    firstLine = true

                # Plot basis sites
                for i in 1:basis_size(unit_cell)
                    basis_pos = unit_cell.basis[:, i]
                    locX = basis_pos[1] 
                    locY = basis_pos[2] 
                    locZ = basis_pos[3] 
                    color = unit_cell.site_colors[i]
                    push!(sitesInfo, [locX, locY, locZ, color])
                end

                # Plot bonds
                for bond in unit_cell.bonds
                    prX = bond.direction[1]
                    prY = bond.direction[2]
                    prZ = bond.direction[3]

                    basis1 = unit_cell.basis[:, bond.site1]
                    basis2 = unit_cell.basis[:, bond.site2]

                    site1X = basis1[1] 
                    site1Y = basis1[2] 
                    site1Z = basis1[3] 

                    site2X = basis2[1] + prX  * pv1[1] + prY  * pv2[1] + prZ  * pv3[1]
                    site2Y = basis2[2] + prX  * pv1[2] + prY  * pv2[2] + prZ  * pv3[2]
                    site2Z = basis2[3] + prX  * pv1[3] + prY  * pv2[3] + prZ  * pv3[3]

                    # Update min/max values
                    maxXValue = max(maxXValue, site1X, site2X)
                    minXValue = min(minXValue, site1X, site2X)
                    maxYValue = max(maxYValue, site1Y, site2Y)
                    minYValue = min(minYValue, site1Y, site2Y)
                    maxZValue = max(maxZValue, site1Z, site2Z)
                    minZValue = min(minZValue, site1Z, site2Z)

                    xData = [site1X, site2X]
                    yData = [site1Y, site2Y]
                    zData = [site1Z, site2Z]
                    bondColor = bond.bond_type 
                    my_palette = colorschemes[:seaborn_bright]

                    plot3d!(unitCellPlot, xData, yData, zData,
                        label="", color=my_palette[bondColor], lw=2)
                end

    # Update bounds based on primitive vectors
    for i in 1:3
        pv = unit_cell.primitive_vectors[:, i]
        minXValue = min(minXValue, pv[1])
        maxXValue = max(maxXValue, pv[1])
        minYValue = min(minYValue, pv[2])
        maxYValue = max(maxYValue, pv[2])
        minZValue = min(minZValue, pv[3])
        maxZValue = max(maxZValue, pv[3])
    end

    # Plot unit cell edges
    for l in unitcellLines
        site1 = l[1]
        site2 = l[2]
        xPoints = [site1[1], site2[1]]
        yPoints = [site1[2], site2[2]]
        zPoints = [site1[3], site2[3]]
        plot3d!(unitCellPlot, xPoints, yPoints, zPoints,
            linestyle=:dash, color="black",
            label=firstLine ? "Unit Cell" : "")
        firstLine = false
    end

    # Set limits and aspect ratio
    plot3d!(unitCellPlot,
        xlims=(minXValue * 1.1, maxXValue * 1.1),
        ylims=(minYValue * 1.1, maxYValue * 1.1),
        zlims=(minZValue * 1.1, maxZValue * 1.1),
        aspect_ratio=:equal,
        xlabel="X", ylabel="Y", zlabel="Z")

    # Plot atoms
    xPositions = getindex.(sitesInfo, 1)
    yPositions = getindex.(sitesInfo, 2)
    zPositions = getindex.(sitesInfo, 3)
    colors = getindex.(sitesInfo, 4)
    paletteAtoms = colorschemes[:Accent_8]    
    for c in unique(colors)
        mask = colors .== c
        colorIndex = Int(c)
            plot3d!(unitCellPlot, xPositions[mask], yPositions[mask], zPositions[mask],
            seriestype=:scatter3d,
            markershape=:circle, 
                   aspect_ratio=:equal,
                   color=paletteAtoms[(colorIndex)],
                   label="Type $colorIndex",
                   markerstrokewidth=2,
                   markersize=4)
    end
    display(unitCellPlot)
    return unitCellPlot
end