"""
    lbp_original(X; kwargs...)
    lbp_original(X, npoints, radius, interpolation=Linear(); kwargs...)
    lbp_original!(out, X, offsets; kwargs...)

Compute the local binary pattern, the original version, of gray image `X` using 3x3
neighborhood matrix.

# Arguments

- `X::AbstractMatrix`: the input image matrix. For colorful images, one can manually convert
  it to some monochrome space, e.g., `Gray` or the L-channel of `Lab`.

[2,4] proposes a generalized interpolation-based version using circular neighborhood matrix,
it produces better result for `rotation=true` case but is usually much slower than the plain
3x3 matrix version. The arguments for this version are:

- `npoints::Int`(4 ≤ npoints ≤ 8): the number of (uniform-spaced) neighborhood points.
- `radius::Real`(radius ≥ 1.0): the radius of the circular.
- `interpolation::Union{Degree, InterpolationType}=Linear()`: the interpolation method used
  to generate non-grid pixel value. In most cases, `Linear()` are good enough. One can also
  try other costly interpolation methods, e.g., `Cubic(Line(OnGrid()))`(also known as
  "bicubic"), `Lanczos()`. See also Interpolations.jl for more choices.

!!! note "neighborhood order differences"
    Different implementation might use different neighborhood orders； this will change the
    encoding result but will not change the overall distribution. For instance,
    `lbp_original(X)` differs from `lbp_original(X, 8, 1, Constant())` only by how `offsets`
    (see below) are ordered; the former uses column-major top-left to bottom-right 3x3 matrix
    order and the latter uses circular order.

Arguments for in-place version:

- `offsets::Tuple`: tuple of neighborhood matrix, e.g., `((0, 1), (0, -1), (1, 0), (-1, 0))`
  specifies the 4-neighborhood matrix. If `X isa Interpolations.AbstractInterpolation`
  holds, then the position values can be float numbers, e.g, `(0.7, 0.7)`.

# Parameters

The parameters control whether and what degree additional encoding passses are used to
compute more patterns that are more robust/invariant to certain changes, e.g., rotation. The
following lists are ordered as encoding order. For example, if `rotation=true` and
`uniform_degree=2`, then rotation encoding will be applied first.

- `rotation=false`: set `true` to generate patterns that are invariant to rotation [3]. For
  example, pattern `0b00001101` is equivalent to `0b01000011` when `rotation=true`.
- `uniform_degree`: the threshold number of pattern uniform degree. From [2] a typical
  choice is `2`.If it is `nothing`(default value) then no uniform encoding is applied.

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

julia> lbp_original(X, 4, 1) # 4-neighbor with circular radius 1
3×3 $(Matrix{UInt8}):
 0x01  0x01  0x00
 0x03  0x02  0x0e
 0x02  0x07  0x00

julia> lbp_original(X; rotation=true)
3×3 $(Matrix{UInt8}):
 0x03  0x01  0x00
 0x0d  0x35  0x1b
 0x05  0x5b  0x00

julia> lbp_original(X; uniform_degree=2)
3×3 $(Matrix{UInt8}):
 0xc0  0x40  0x00
 0x09  0x09  0x09
 0x09  0x09  0x00
```

# Extended help

## Local binary pattern

The following is how local binary pattern is calculated, the original version[1]:

```text
3x3 block     center-thresholded     weights         multiplied by weights      sum
6  7  9         1  1  1              1  8  32          1   8  32
5  6  3  ==>    0  x  0     ==>      2  x  64  ==>     0   x  0            ==>  169
2  1  7         1  0  1              4  16 128         0   0  128
```

Any binary pattern of length 8, i.e., the center-thresholded result, can be uniquely
represented as an `UInt8` value; the weighted sum is the encoding process.

## Rotation-invariant encoding

The rotation-invariant encoding is to map all elements in the bitrotation equivalent class
to the minimal value of this class. For example, `0b11010000` and `0b01000011` belongs to
the same class because `bitrotate(0b01000011, -2) == 0b11010000`, thus both values are
mapped to `0b00001101`. See also Eq.(8) in [2].

For 3x3 neighborhood matrix, applying rotation-invariant encoding decreases the possible
number of binary patterns from ``256`` to ``36``.

The interpolation-based version provides more robust result for rotation-invariant pattern,
see [2,4] for more details.

## Uniform encoding

Authors of [2] states that certain local binary patterns are fundamental properties of
texture, providing the vast majority, sometimes over 90 percent, of all 3x3 patterns. Those
patterns are called "uniform" as they contain very few spatial transitions. Uniform degree
is an additional encoding pass that controls at what circumstances can we set the block to
miscellaneous class.

For example, if `uniform_degree=2`, then `0b00001101` will be encoded as `9` (type
miscellaneous) because it has ``3`` bit transitions, and `0b00001100` will be unchanged
because it only has ``2`` bit transitions.

# References

- [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
- [2] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
- [3] Pietikäinen, Matti, Timo Ojala, and Zelin Xu. "Rotation-invariant texture classification using feature distributions." _Pattern recognition_ 33.1 (2000): 43-52.
- [4] T. Ojala, M. Pietikainen, and T. Maenpaa, “Multiresolution gray-scale and rotation invariant texture classification with local binary patterns,” _IEEE Trans. Pattern Anal. Machine Intell._, vol. 24, no. 7, pp. 971–987, Jul. 2002, doi: 10.1109/TPAMI.2002.1017623.
"""
function lbp_original(X::AbstractArray; kwargs...)
    # The original version [1] uses clockwise order; here we use anti-clockwise order
    # because Julia is column-major order. If we consider memory order differences then they
    # are equivalent.
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    lbp_original!(similar(X, UInt8), X, offsets; kwargs...)
end
function lbp_original!(
        out,
        X::AbstractMatrix{T},
        offsets::Tuple;
        rotation::Bool=false,
        uniform_degree::Union{Nothing,Int}=nothing,
        ) where T<:Union{Real,Gray}
    length(offsets) > 8 && throw(ArgumentError("length(offsets) >= 8 is not supported."))
    outerR = CartesianIndices(X)
    r = CartesianIndex(
        ceil(Int, maximum(abs.(extrema(first.(offsets))))),
        ceil(Int, maximum(abs.(extrema(last.(offsets)))))
    )
    innerR = first(outerR)+r:last(outerR)-r

    # TODO(johnnychen94): use LoopVectorization
    @inbounds for I in innerR
        gc = X[I]
        # This inner loop fuses the binary pattern build stage and encoding stage for
        # better performance.
        rst = 0
        for i in 1:length(offsets)
            p = I.I .+ offsets[i]
            rst += ifelse(gc <= _inbounds_getindex(X, p), 1, 0) << (i-1)
        end
        out[I] = rst
    end

    # boundary conditions are undefined in [1]; here we directly skip out-of-boundary values
    @inbounds for I in EdgeIterator(outerR, innerR)
        gc = X[I]
        rst = 0
        for i in 1:length(offsets)
            p = I.I .+ offsets[i]
            checkbounds(Bool, X, p...) || continue
            rst += ifelse(gc <= _inbounds_getindex(X, p), 1, 0) << (i-1)
        end
        out[I] = rst
    end

    # The building is cached and chained(if there are multiple encoding passes) thus the
    # cost is decreased to one brocasting to `getindex` and `setindex!`.
    encoding_table = build_LBP_encoding_table(UInt8; rotation=rotation, uniform_degree=uniform_degree)
    if !isnothing(encoding_table)
        @. out = encoding_table[out + 1]
    end

    return out
end

function lbp_original(X::AbstractArray, npoints, radius, interpolation=Linear(); kwargs...)
    # TODO(johnnychen94): support npoints=24 as [2, 4] indicate
    4<= npoints <= 8 || throw(ArgumentError("currently only support 4<= npoints <=8, instead it is $(npoints)."))
    interpolation = wrap_BSpline(interpolation)
    offsets = _circular_neighbor_offsets(npoints, radius)
    if interpolation == BSpline(Constant())
        offsets = map(offsets) do o
            round.(Int, o, RoundToZero)
        end
        itp = X
    else
        itp = interpolate(X, interpolation)
    end
    lbp_original!(similar(X, UInt8), itp, offsets; kwargs...)
end
