struct ExpansionVertices{V<:Integer} <: AbstractVertices{V}
    vertices::AbstractSet{V}
end

function ExpansionVertices()
    ExpansionVertices(BitSet())
end

function ExpansionVertices(i::Int)
    ExpansionVertices(BitSet(i))
end

function ExpansionVertices{V}(i::Int) where {V<:Integer}
    ExpansionVertices{V}(BitSet(i))
end

function ExpansionVertices(itr::AbstractVector{<:Integer})
    ExpansionVertices(BitSet(itr))
end

vertices(vs::ExpansionVertices) = vs.vertices
Base.sort(vs::ExpansionVertices) = vs
Base.collect(vs::ExpansionVertices) = collect(vs.vertices)
Base.intersect(vs1::ExpansionVertices, vs2::ExpansionVertices) = ExpansionVertices(intersect(vertices(vs1), vertices(vs2)))
Base.setdiff(vs1::ExpansionVertices, vs2::ExpansionVertices) = ExpansionVertices(setdiff(vertices(vs1), vertices(vs2)))
Base.union(vs1::ExpansionVertices, vs2::ExpansionVertices) = ExpansionVertices(union(vertices(vs1), vertices(vs2)))
Base.in(v::V, vs::ExpansionVertices{V}) where {V<:Integer} = v in vertices(vs)
Base.eltype(vs::ExpansionVertices{V}) where {V} = V
