module Vertices

import LINCEGE:
    _NI

abstract type AbstractVertices{V} end

vertices(vs::AbstractVertices) = _NI("vertices")
Base.collect(vs::AbstractVertices) = _NI("Base.collect")
Base.sort(vs::AbstractVertices) = _NI("Base.sort")
Base.intersect(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.intersect")
Base.setdiff(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.setdiff")
Base.union(vs1::AbstractVertices, vs2::AbstractVertices) = _NI("Base.union")
Base.in(v, vs::AbstractVertices) = _NI("Base.in")
Base.eltype(vs::AbstractVertices) = _NI("Base.eltype")

Base.length(vs::AbstractVertices) = length(vertices(vs))

Base.isequal(vs1::AbstractVertices, vs2::AbstractVertices) = vs1 == vs2
Base.:(==)(vs1::AbstractVertices, vs2::AbstractVertices) = (collect(vs1) == collect(vs2))
Base.hash(vs::AbstractVertices, h::UInt) = hash(collect(sort(vs)), h)
Base.contains(vs::AbstractVertices, v) = v in vs
Base.haskey(vs::AbstractVertices, v) = v in vs

function Base.union(vs::AbstractVertices, itr)
    vstemp = vs
    for x in itr
        vstemp = union(vstemp, x)
    end
    vstemp
end

Base.getindex(vec::Vector, vs::AbstractVertices) = vec[collect(vs)]
Base.getindex(mat::Matrix, vs1::AbstractVertices, vs2::AbstractVertices) = mat[collect(vs1), collect(vs2)]
Base.getindex(mat::Matrix, vs1::AbstractVertices, i::Int) = mat[collect(vs1), i]
Base.getindex(mat::Matrix, i::Int, vs1::AbstractVertices) = mat[i, collect(vs1)]

Base.iterate(vs::AbstractVertices) = iterate(vertices(vs))
Base.iterate(vs::AbstractVertices, state) = iterate(vertices(vs), state)

Base.show(io::IO, vs::AbstractVertices) = print(io, "Vertices: ", collect(vs))

include("ExpansionVertices.jl")
include("LatticeVertices.jl")

export AbstractVertices,
    ExpansionVertices,
    LatticeVertices
end
