module UnitCells
import LINCEGE:
    _NI

abstract type AbstractUnitCell end

basis_size(unit_cell::AbstractUnitCell) = _NI("basis_size")
dimension(unit_cell::AbstractUnitCell) = _NI("dimension")
shift_unit_cell(unit_cell::AbstractUnitCell, shift_vector::AbstractVector{<:Int}) = _NI("shift_unit_cell")
# matrix input is of the form [shift_x shift_y shift_z; ...] where each row is a shift vector for the unit cell
shift_unit_cell(unit_cell::AbstractUnitCell, shift_vectors::AbstractMatrix{<:Int}) = _NI("shift_unit_cell")
image_unit_cell(unit_cell::AbstractUnitCell) = _NI("image_unit_cell")

include("UnitCell.jl")

end
