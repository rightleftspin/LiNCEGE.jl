struct TranslationHash{H<:Unsigned} <: AbstractGraphHash{H}
    hash::H
end

function TranslationHash(vs::AbstractVertices, lattice::AbstractLattice)
    TranslationHash(hash(translation_form(lattice, vs)))
end
