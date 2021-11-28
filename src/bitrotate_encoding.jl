const _RotationEncodingTables = Dict()

"""
    rotation_encoding_table(T, nbits)
    rotation_encoding_table(T, X::AbstractUnitRange)

Build the rotation-invariant encoding table. If global runtime cache exists, then reuse it.
"""
rotation_encoding_table(::Type{T}, nbits::Int) where T = rotation_encoding_table(T, 0:2^nbits-1)
function rotation_encoding_table(::Type{T}, X::AbstractUnitRange) where T
    d = _RotationEncodingTables
    k = (T, minimum(X), maximum(X))
    haskey(d, k) && return d[k]
    lookup = d[k] = RotationEncodingTable{T}(X)
    return lookup
end

"""
    RotationEncodingTable{T<:Unsigned}(X)

The lookup table for the quotient space constructed using `bitrotate` equivalance relation
on bit representation of datatype `T`.

By definition, ``a`` is equivalant to ``b``(``a ~ b``) if there exists ``n`` such that
`bitrotate(a, n) == b`.

This equivalance class gives a unique encoding map from ``x∈X`` to ``a``, where ``a`` is the
minimal value of the equivalance class ``[x]``. For instance, the equivalance class of
`0b10101001` (UInt8) is

```julia
0b10101001 # 169
0b01010011 # 83
0b10100110 # 166
0b01001101 # 77
0b10011010 # 154
0b00110101 # 53
0b01101010 # 106
0b11010100 # 212
```

thus `f(169) == 0x35 == 53`.

```jldoctest
julia> using LocalBinaryPatterns: RotationEncodingTable

julia> X = RotationEncodingTable{UInt8}(0:255);

julia> X[169] # 53
0x35

julia> X = RotationEncodingTable{UInt8}(128:255);

julia> X[169] # 154
0x9a
```

The values in the lookup table is implemented in lazy cache manner so the real computation
only happens once at the first retrieval.

This lookup table is used to build an efficient implementation of rotation-invariant local
binary pattern [1].

# References

- [1] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
"""
struct RotationEncodingTable{T<:Unsigned,AT<:AbstractUnitRange,LAT} <: AbstractVector{T}
    X::AT
    Y::LAT
end
function RotationEncodingTable{T}(X::AbstractUnitRange) where T
    if !((maximum(X) <= typemax(T)) && (minimum(X) >= typemin(T)))
        throw(ArgumentError("range $(first(X)):$(last(X)) can't be uniquely represented by type $T."))
    end

    if log2(length(X)) >= 26
        datasize = round(length(X)/1024/1024/1024, digits=3)
        @warn "$(datasize)GB data will be created."
    end

    Y = Vector{Union{T,Missing}}(missing, length(X))
    RotationEncodingTable{T,typeof(X),typeof(Y)}(X, Y)
end

@inline Base.axes(X::RotationEncodingTable) = (first(X.X):last(X.X), )
@inline Base.size(X::RotationEncodingTable) = (length(X.X)-1)
@inline function Base.getindex(X::RotationEncodingTable{T}, x::Int) where T
    i = x-first(axes(X, 1))+1
    v = X.Y[i]
    if ismissing(v)
        v = minimum_bitrotate(T(x), X.X)
        X.Y[i] = v
    end
    return T(v)
end

function minimum_bitrotate(x::T, X::AbstractUnitRange) where T<:Unsigned
    function bit_rotate_project_to(x, n, X)
        v = bitrotate(x, n)
        return v in X ? v : x
    end
    minimum(bit_rotate_project_to(x, n, X) for n in 0:8sizeof(T)-1)
end
