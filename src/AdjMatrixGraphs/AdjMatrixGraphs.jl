abstract type AbstractAdjMatrixGraph end

nv(g::AbstractAdjMatrixGraph) = size(g.adj_matrix, 1)
nw(g::AbstractAdjMatrixGraph) = length(g.weight_info)

labels(g::AbstractAdjMatrixGraph) = diag(g.adj_matrix)
adj_matrix(g::AbstractAdjMatrixGraph) = tril(g.adj_matrix, -1) + tril(g.adj_matrix, 1)

edge_list(g::AbstractAdjMatrixGraph) = adj_matrix_to_edge_list(adj_matrix(g))


Base.show(io::IO, g::AbstractAdjMatrixGraph) = print(
    io,
    "Graph with $(nv(g)) vertices, $(length(edge_list(g))) bonds, and $(nw(g)) unique weights",
)
