const _UniformEncodingTables = Dict()
"""
    uniform_encoding_table(T, X::AbstractUnitRange, degree::Int)

Get the encoding table of `X` using uniform degree filter.

A uniform degree filter is defined as `f(x) = U(x) > degree ? P+1 : x`, where `P` is number
of bits of `x`, and `U(x)` is the circular bit transition of `x` with bit representation of
datatype `T`.

For example, `0b0000001` has ``2`` bit transitions thus is unchanged, while `0b00000101` has
``4`` bit transitions. The transition is count in circular sense, i.e., the highest and
lowest bits are also compared.

```jldoctest
julia> using LocalBinaryPatterns: uniform_encoding_table

julia> X = uniform_encoding_table(UInt8, 0:255, 2);

julia> X[0b0000001]
0x01

julia> X[0b0000101] # miscellaneous
0x09
```

This function is used to distinguish local binary patterns from texture patterns and
miscellaneous patterns. See [1] for more information.

!!! info "Runtime caches"
    For better performance, this function uses a runtime cache to store the lookup and
    shared among all callers. The result is expected to be used in read-only mode.

# References

- [1] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
"""
function uniform_encoding_table(::Type{T}, X::AbstractUnitRange, degree::Int) where T
    d = _UniformEncodingTables
    k = (T, degree, minimum(X), maximum(X))
    haskey(d, k) && return d[k]
    lookup = d[k] = _uniform_encoding_table(T, X, degree)
    return lookup
end
uniform_encoding_table(::Type{T}, nbits::Int, degree::Int) where T = uniform_encoding_table(T, 0:2^nbits-1, degree)

function _uniform_encoding_table(::Type{T}, X::AbstractUnitRange, degree::Int) where T<:Unsigned
    if !((maximum(X) <= typemax(T)) && (minimum(X) >= typemin(T)))
        throw(ArgumentError("range $(first(X)):$(last(X)) can't be uniquely represented by type $T."))
    end

    function _count_bit_transitions(x)
        # Eq. (10) in Ojala 2002
        count_ones(x << 7 & typemax(typeof(x))) + count_ones(x ⊻ (x >> 1))
    end
    # Eq. (9) for Ojala 2002
    lookup = map(X) do x
        n = _count_bit_transitions(T(x))
        T(ifelse(n > degree, 8sizeof(T)+1, x))
    end
    OffsetVector(lookup, X)
end
