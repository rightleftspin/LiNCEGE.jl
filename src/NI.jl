"""
    NotImplementedError{M}(m)

`Exception` thrown when a method required for a NLCE is not implemented.
Taken directly from Graphs.jl.
"""
struct NotImplementedError{M} <: Exception
    m::M
    NotImplementedError(m::M) where {M} = new{M}(m)
end

function Base.showerror(io::IO, ie::NotImplementedError)
    return print(io, "method $(ie.m) not implemented.")
end

_NI(m) = throw(NotImplementedError(m))
