# This implements the solver for Eq. (8) in the Ojala 2002 as lookup table.
# Computing the encoding table using naive implementation is time-consuming,
# since it's read-only, we simply cache the encoding table.
const _rotation_invariant_encoding_tables = Dict()
const SupportedEncodingTypes = Union{UInt8, UInt16}
function build_rotation_invariant_encoding_table(::Type{T}) where T<:SupportedEncodingTypes
    d = _rotation_invariant_encoding_tables
    haskey(d, T) && return d[T]
    d[T] = _freeze(T, _build_inverse_table(_bitrotate_quotation_space(T)))
    return d[T]
end

# Mathematically, this is the quotation space under circular bitshift of the N-bits binary
# pattern space. For instance, 0b00001101 and 0b01000011 belong to the same equivalance
# class.
# TODO(johnnychen94): maybe support UInt32 and beyond by providing a more efficient implementation
function _bitrotate_quotation_space(::Type{T}) where T<:SupportedEncodingTypes
    s = Vector{T}[]
    for x in typemin(T):typemax(T)
        # without the following skipping mechanism actually runs faster
        # any(c->x in c, s) && continue
        push!(s, sort!(bitrotate.(x, 0:8sizeof(T)-1)))
    end
    Dict(minimum(c) => c for c in unique!(s))
end
function _build_inverse_table(d::Dict{T,<:AbstractVector{T}}) where {T}
    id = Vector{T}(undef, typemax(T)-typemin(T)+1)
    for k in keys(d)
        for q in d[k]
            id[q+1] = k
        end
    end
    return id
end
_freeze(::Type{T}, v::Vector) where T = SVector{typemax(T)-typemin(T)+1}(v)
