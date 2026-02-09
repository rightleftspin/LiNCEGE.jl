function to_matrix(vectors::AbstractVector{<:AbstractVector{<:T}}) where {T}
    num_rows = length(vectors)
    num_cols = length(vectors[1])
    matrix = Matrix{T}(undef, num_rows, num_cols)
    for i in 1:num_rows
        if length(vectors[i]) != num_cols
            throw(DimensionMismatch("All vectors must have the same length"))
        end
        matrix[i, :] = vectors[i]
    end

    matrix
end
