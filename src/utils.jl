const SupportedEncodingTypes = Union{UInt8, UInt16}
const _LBP_encoding_table = Dict()
function build_LBP_encoding_table(::Type{T};
        rotation::Bool,
        uniform_degree::Union{Nothing,Int}=nothing
    ) where T<:SupportedEncodingTypes

    d = _LBP_encoding_table
    p = (T, rotation, uniform_degree)
    haskey(d, p) && return d[p]

    identity_lookup = identity.(typemin(T):typemax(T))
    rot_lookup = rotation ? _build_inverse_table(_bitrotate_quotation_space(T)) : identity_lookup
    uniform_lookup = !isnothing(uniform_degree) ? _uniform_encoding_table(T, uniform_degree) : identity_lookup

    # Chaining multiple encoding passes into one lookup table so that we can move as
    # much computation to warm-up phase as we can.
    lookup = d[p] = _freeze(T, uniform_lookup[rot_lookup.+1])
    return lookup == identity_lookup ? nothing : lookup
end

function _bitrotate_quotation_space(::Type{T}) where T<:SupportedEncodingTypes
    # Mathematically, this is the quotation space under circular bitshift of the N-bits binary
    # pattern space. For instance, 0b00001101 and 0b01000011 belong to the same equivalance
    # class.
    # TODO(johnnychen94): maybe support UInt32 and beyond by providing a more efficient implementation

    # This implements the solver for Eq. (8) in the Ojala 2002 as lookup table.
    # Computing the encoding table using naive implementation is time-consuming,
    # since it's read-only, we simply cache the encoding table.
    s = Vector{T}[]
    for x in typemin(T):typemax(T)
        # without the following skipping mechanism actually runs faster
        # any(c->x in c, s) && continue
        push!(s, sort!(bitrotate.(x, 0:8sizeof(T)-1)))
    end
    Dict(minimum(c) => c for c in unique!(s))
end

function _uniform_encoding_table(::Type{T}, degree::Int) where T<:SupportedEncodingTypes
    function _count_bit_transitions(x)
        count_ones(x << 7 & typemax(typeof(x))) + count_ones(x ⊻ (x >> 1))
    end
    # Eq. (9) for Ojala 2002
    map(typemin(T):typemax(T)) do x
        n = _count_bit_transitions(x)
        ifelse(n > degree, 8sizeof(T)+1, x)
    end
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
