module UnitCells

import LINCEGE:
    _NI

abstract type AbstractUnitCell end

shift_unit_cell(unit_cell::AbstractUnitCell, shift_vector::Vector{Int}) = _NI("shift_unit_cell")
image_unit_cell(unit_cell::AbstractUnitCell) = _NI("image_unit_cell")

end
