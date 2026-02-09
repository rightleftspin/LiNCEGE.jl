struct ClusterExpansion{H<:AbstractGraphHash,C<:AbstractCluster} <: AbstractClusterExpansion{H,C}
    clusters::Dict{H,C}
    multiplicities::Dict{H,AbstractVector{Float64}}
    subgraphs::Dict{H,Subgraphs}
end

function SiteExpansion(clusters::IsomorphicClusters, lattice::AbstractClusterExpansionLattice)
    maxOrder = max_order(lattice)
    clustersDict = Dict{IsomorphicHash,AbstractCluster}()
    multiplicitiesDict = Dict{IsomorphicHash,Vector{<:Real}}()
    subgraphsDict = Dict{IsomorphicHash,Subgraphs}()

    @info "Starting subgraph enumeration for $(length(clusters)) clusters"
    for (ghash, cluster) in clusters
        clustersDict[ghash] = cluster
        multiplicitiesDict[ghash] = zeros(Float64, maxOrder + 1)
        subgraphsDict[ghash] = Subgraphs(cluster, lattice)
    end
    @info "Finished subgraph enumeration"

    SiteExpansion{IsomorphicHash,IsomorphicCluster}(clustersDict, multiplicitiesDict, subgraphsDict)
end

Base.getindex(ce::SiteExpansion, cluster::AbstractCluster, order::Int) = @inbounds ce.multiplicities[ghash(cluster)][order]
Base.getindex(ce::SiteExpansion, ghash::AbstractGraphHash, order::Int) = @inbounds ce.multiplicities[ghash][order]
Base.setindex!(ce::SiteExpansion, lattice_constant::Real, cluster::AbstractCluster, order::Int) = @inbounds ce.multiplicities[ghash(cluster)][order] = lattice_constant
Base.setindex!(ce::SiteExpansion, lattice_constant::Real, ghash::AbstractGraphHash, order::Int) = @inbounds ce.multiplicities[ghash][order] = lattice_constant
subgraphs(ce::SiteExpansion, cluster::AbstractCluster) = @inbounds ce.subgraphs[ghash(cluster)]
subgraphs(ce::SiteExpansion, ghash::AbstractGraphHash) = @inbounds ce.subgraphs[ghash]
Base.iterate(ce::SiteExpansion) = iterate(zip(keys(ce.clusters), values(ce.clusters), values(ce.multiplicities)))
Base.iterate(ce::SiteExpansion, state) = iterate(zip(keys(ce.clusters), values(ce.clusters), values(ce.multiplicities)), state)
