"""
    lbp_original(X)
    lbp_original!(out, X)

Compute the local binary pattern, the original version, of gray image using 3x3
neighborhood matrix.

```text
3x3 block     center-thresholded     weights         multiplied by weights      sum
6  7  9         1  1  1              1  8  32          1   8  32
5  6  3  ==>    0  x  0     ==>      2  x  64  ==>     0   x  0            ==>  169
2  1  7         1  0  1              4  16 128         0   0  128
```

# Examples

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> X = [6 7 9; 5 6 3; 2 1 7]
3×3 $(Matrix{Int}):
 6  7  9
 5  6  3
 2  1  7

julia> lbp_original(X)
3×3 $(Matrix{UInt8}):
 0xc0  0x40  0x00
 0x68  0xa9  0x1b
 0x28  0x6b  0x00
```

# References

- [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
- [2] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
"""
lbp_original(X::AbstractArray) = lbp_original!(similar(X, UInt8), X)
function lbp_original!(out, X::AbstractMatrix{T}) where T<:Union{Real,Gray}
    # nearest interpolation, 3x3 neighborhood

    # The original version [1] uses clockwise order; here we use anti-clockwise order
    # because Julia is column-major order. If we consider memory order differences then they
    # are equivalent.
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    offsets = CartesianIndex.(offsets)
    outerR = CartesianIndices(X)
    innerR = first(outerR)+oneunit(first(outerR)):last(outerR)-oneunit(last(outerR))

    # TODO(johnnychen94): use LoopVectorization
    @inbounds for I in innerR
        gc = X[I]
        # This inner loop fuses the binary pattern build stage and encoding stage for
        # better performance.
        rst = 0
        for i in 1:length(offsets)
            o = offsets[i]
            rst += ifelse(gc <= X[I+o], 1, 0) << (i-1)
        end
        out[I] = rst
    end

    # boundary conditions are undefined in [1]; here we directly skip out-of-boundary values
    @inbounds for I in EdgeIterator(outerR, innerR)
        gc = X[I]
        rst = 0
        for i in 1:length(offsets)
            o = offsets[i]
            checkbounds(Bool, X, I+o) || continue
            rst += ifelse(gc <= X[I+o], 1, 0) << (i-1)
        end
        out[I] = rst
    end

    return out
end
