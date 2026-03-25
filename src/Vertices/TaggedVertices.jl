struct TaggedVertices{Tag,V<:Integer} <: AbstractVertices{V}
        bitset::BitSet
end

vertices(vs::TaggedVertices) = vs.bitset

Base.sort(vs::TaggedVertices) = vs
Base.collect(vs::TaggedVertices) = collect(vs.bitset)

Base.intersect(vs1::T, vs2::T) where {T<:TaggedVertices} = T(intersect(vs1.bitset, vs2.bitset))
Base.setdiff(vs1::T, vs2::T) where {T<:TaggedVertices} = T(setdiff(vs1.bitset, vs2.bitset))
Base.union(vs1::T, vs2::T) where {T<:TaggedVertices} = T(union(vs1.bitset, vs2.bitset))

Base.in(v::Integer, vs::TaggedVertices) = v in vs.bitset
Base.eltype(::TaggedVertices{Tag,V}) where {Tag,V} = V

TaggedVertices{Tag,V}() where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet())
TaggedVertices{Tag,V}(i::Int) where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet(i))
TaggedVertices{Tag,V}(itr::AbstractVector{<:Integer}) where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet(itr))

const ExpansionVertices{V} = TaggedVertices{:expansion,V}

ExpansionVertices() = ExpansionVertices{Int}(BitSet())
ExpansionVertices(i::Int) = ExpansionVertices{Int}(BitSet(i))
ExpansionVertices(itr::AbstractVector{<:Integer}) = ExpansionVertices{Int}(BitSet(itr))
ExpansionVertices(bs::BitSet) = ExpansionVertices{Int}(bs)

const LatticeVertices{V} = TaggedVertices{:lattice,V}

LatticeVertices() = LatticeVertices{Int}(BitSet())
LatticeVertices(i::Int) = LatticeVertices{Int}(BitSet(i))
LatticeVertices(itr::AbstractVector{<:Integer}) = LatticeVertices{Int}(BitSet(itr))
LatticeVertices(bs::BitSet) = LatticeVertices{Int}(bs)
