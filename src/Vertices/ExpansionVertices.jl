struct ExpansionVertices{V<:AbstractSet{<:Unsigned}}
        vertices::V
end

ExpansionVertices(v::Vertex) where {Vertex<:Unsigned} = ExpansionVertices(BitSet(v))

Base.union(exp_v1::ExpansionVertices, v::Vertex) where {Vertex<:Unsigned} = ExpansionVertices(union(exp_v1.vertices, v))

Base.union(exp_v1::ExpansionVertices, exp_v2::ExpansionVertices) = ExpansionVertices(union(exp_v1.vertices, exp_v2.vertices))

function Base.union(exp_v1::ExpansionVertices, itr)
        exp_temp = exp_v1
        for x in itr
                exp_temp = union(exp_temp, x)
        end

        exp_temp
end

Base.intersect(exp_v1::ExpansionVertices, exp_v2::ExpansionVertices) = ExpansionVertices(intersect(exp_v1, exp_v2))

Base.collect(exp_v::ExpansionVertices) = collect(exp_v.vertices)
Base.in(v, exp_v::ExpansionVertices) = v in exp_v.vertices

Base.length(exp_v::ExpansionVertices) = length(exp_v.vertices)

Base.getindex(vec::AbstractVector, exp_v::ExpansionVertices) = vec[collect(exp_v)]

Base.pop!(exp_v::ExpansionVertices) = pop!(exp_v.vertices)




