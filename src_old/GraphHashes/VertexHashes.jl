struct VertexHash{H<:Unsigned} <: AbstractGraphHash{H}
    hash::H
end

VertexHash(vs::AbstractVertices) = VertexHash(hash(vs))
