
struct BondCluster{V<:ExpansionVertices,H<:Unsigned} <: AbstractCluster
        expansion_vertices::V
        ghash::H
end
