
struct DistanceCluster{V<:ExpansionVertices,H<:Unsigned} <: AbstractCluster
    expansion_vertices::V
    ghash::H
end
