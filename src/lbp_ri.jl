"""
    lbp_ri(X)
    lbp_ri!(out, X)

Compute the local binary pattern, the rotation invariant version, of gray image using 3x3
neighborhood matrix.

```text
3x3 block     center-thresholded     weights         multiplied by weights      sum      rotation invariant encoding
6  7  9         1  1  1              1  8  32          1   8  32
5  6  3  ==>    0  x  0     ==>      2  x  64  ==>     0   x  0            ==>  169  ==>    53
2  1  7         1  0  1              4  16 128         0   0  128
```

The rotation encoding is to map all elements in the bitrotation equivalent class to the
minimal value of this class. For example, `0b11010000` and `0b01000011` belongs to the same
class because `bitrotate(0b01000011, -2) == 0b11010000`, thus both values are mapped to
`0b00001101`, which is `minimum(bitrotate.(0b11010000, 0:7))`.

# Examples

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> X = [6 7 9; 5 6 3; 2 1 7]
3×3 $(Matrix{Int}):
 6  7  9
 5  6  3
 2  1  7

julia> lbp_ri(X)
3×3 $(Matrix{UInt8}):
 0x03  0x01  0x00
 0x0d  0x35  0x1b
 0x05  0x5b  0x00
```

# References

- [1] Pietikäinen, Matti, Timo Ojala, and Zelin Xu. "Rotation-invariant texture classification using feature distributions." _Pattern recognition_ 33.1 (2000): 43-52.
- [2] T. Ojala, M. Pietikainen, and T. Maenpaa, “Multiresolution gray-scale and rotation invariant texture classification with local binary patterns,” _IEEE Trans. Pattern Anal. Machine Intell._, vol. 24, no. 7, pp. 971–987, Jul. 2002, doi: 10.1109/TPAMI.2002.1017623.
"""
lbp_ri(X::AbstractArray) = lbp_ri!(similar(X, UInt8), X)

function lbp_ri!(out, X::AbstractMatrix{T}) where T<:Union{Real, Gray}
    # nearest interpolation, 3x3 neighborhood, rotation invariant
    encoding_table = build_circular_invariant_encoding_table(UInt8)
    lbp_original!(out, X)
    @. out = encoding_table[out + 1]
    return out
end
