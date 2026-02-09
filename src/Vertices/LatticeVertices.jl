struct LatticeVertices{V<:Integer} <: AbstractVertices{V}
    vertices::AbstractSet{V}
end

function LatticeVertices()
    LatticeVertices(BitSet())
end

function LatticeVertices(i::Int)
    LatticeVertices(BitSet(i))
end

function LatticeVertices{V}(i::Int) where {V<:Integer}
    LatticeVertices{V}(BitSet(i))
end

function LatticeVertices(itr::AbstractVector{<:Integer})
    LatticeVertices(BitSet(itr))
end

vertices(vs::LatticeVertices) = vs.vertices
Base.sort(vs::LatticeVertices) = vs
Base.collect(vs::LatticeVertices) = collect(vs.vertices)
Base.intersect(vs1::LatticeVertices, vs2::LatticeVertices) = LatticeVertices(intersect(vertices(vs1), vertices(vs2)))
Base.setdiff(vs1::LatticeVertices, vs2::LatticeVertices) = LatticeVertices(setdiff(vertices(vs1), vertices(vs2)))
Base.union(vs1::LatticeVertices, vs2::LatticeVertices) = LatticeVertices(union(vertices(vs1), vertices(vs2)))
Base.in(v::V, vs::LatticeVertices{V}) where {V<:Integer} = v in vertices(vs)
Base.eltype(vs::LatticeVertices{V}) where {V} = V
