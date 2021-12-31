###
# Various LBP modes
###

abstract type LocalBinaryPattern end

struct LBP{N,F} <: LocalBinaryPattern
    mode::F # mode(X::AbstractArray{T,N}, I::CartesianIndex{N}, offsets::Tuple)
    block_size::Dims{N}
    function LBP(mode::F, block_size::Dims{N}) where {N,F}
        all(block_size .> 0) || throw(ArgumentError("block size `block_size=$(block_size)` should be positive integers."))
        new{N, F}(mode, block_size)
    end
end
LBP{N}(mode, block_size::Int=1) where N = LBP(mode, ntuple(_->block_size, N))

Base.@propagate_inbounds _default_lbp_mode(X, I, offsets) = X[I]
Base.@propagate_inbounds function _default_lbp_mode(X::IntegralArray, I, bo, offsets)
    Rp = I:min(I+bo, last(CartesianIndices(X)))
    X[first(Rp)..last(Rp)]/length(Rp)
end

Base.@propagate_inbounds function lbp_encode(lbp::LBP, X::AbstractMatrix, I::CartesianIndex{2}, offsets::Tuple)
    rst = 0
    gc = lbp.mode(X, I, offsets)
    for i in 1:length(offsets)
        p = I.I .+ offsets[i]
        @boundscheck checkbounds(Bool, X, p...) || continue
        rst += ifelse(gc <= _inbounds_getindex(X, p), 1, 0) << (i-1)
    end
    return rst
end


Base.@propagate_inbounds function mlbp_encode_inbounds(
        lbp::LBP,
        iX::IntegralArray,
        I::CartesianIndex{2},
        bo::CartesianIndex{2},
        offsets::Tuple)
    rst = 0
    gc = lbp.mode(iX, I, bo, offsets)
    n = prod(bo.I .+ 1)
    for i in 1:length(offsets)
        p = CartesianIndex(I.I .+ offsets[i])
        rst += ifelse(gc <= iX[p..p+bo]/n, 1, 0) << (i-1)
    end
    return rst
end
function mlbp_encode(
        lbp::LBP,
        iX::IntegralArray,
        I::CartesianIndex{2},
        bo::CartesianIndex{2},
        offsets::Tuple)
    rst = 0
    Rlast = last(CartesianIndices(iX))
    gc = lbp.mode(iX, I, bo, offsets)
    for i in 1:length(offsets)
        p = CartesianIndex(I.I .+ offsets[i])
        checkbounds(Bool, iX, p.I...) || continue
        Rp = p:min(p+bo, Rlast)
        gp = iX[first(Rp)..last(Rp)]/length(Rp)
        rst += ifelse(gc <= gp, 1, 0) << (i-1)
    end
    return rst
end


# non-interpolation version
"""
    local_binary_pattern([f], X; [rotation], [uniform_degree])

Compute the local binary pattern of gray image `X`.

# Arguments

- `f`: function `f(X, I, offsets)` computes the block mode value. In the original version
  [1] it directly uses the center pixel value, i.e., `f(X, I, offsets) = X[I]`.
  See also [`average_mode`](@ref LocalBinaryPatterns).
- `X::AbstractMatrix`: the input image matrix. For colorful images, one can manually convert
  it to some monochrome space, e.g., `Gray`, the L-channel of `Lab`. One could also do
  channelwise LBP and then concatenate together.

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

julia> local_binary_pattern(X)
3×3 $(Matrix{UInt8}):
 0xc0  0x40  0x00
 0x68  0xa9  0x1b
 0x28  0x6b  0x00

julia> local_binary_pattern(X; rotation=true)
3×3 $(Matrix{UInt8}):
 0x03  0x01  0x00
 0x0d  0x35  0x1b
 0x05  0x5b  0x00

julia> local_binary_pattern(X; uniform_degree=2)
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

For example, in 8-bit mode, if `uniform_degree=2`, then `0b00001101` will be encoded as `9`
(type miscellaneous) because it has ``3`` bit transitions, and `0b00001100` will be
unchanged because it only has ``2`` bit transitions.

# References

- [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
- [2] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
- [3] Pietikäinen, Matti, Timo Ojala, and Zelin Xu. "Rotation-invariant texture classification using feature distributions." _Pattern recognition_ 33.1 (2000): 43-52.
- [4] T. Ojala, M. Pietikainen, and T. Maenpaa, “Multiresolution gray-scale and rotation invariant texture classification with local binary patterns,” _IEEE Trans. Pattern Anal. Machine Intell._, vol. 24, no. 7, pp. 971–987, Jul. 2002, doi: 10.1109/TPAMI.2002.1017623.
"""
function local_binary_pattern(f, X::AbstractArray; rotation::Bool=false, uniform_degree::Union{Nothing,Int}=nothing)
    # The original version [1] uses clockwise order; here we use anti-clockwise order
    # because Julia is column-major order. If we consider memory order differences then they
    # are equivalent.
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    lookups = _build_local_binary_pattern_lookup(UInt8, 8; rotation=rotation, uniform_degree=uniform_degree)
    local_binary_pattern!(similar(X, UInt8), LBP{ndims(X)}(f), X, offsets, lookups)
end
local_binary_pattern(X::AbstractArray; kwargs...) = local_binary_pattern(_default_lbp_mode, X; kwargs...)

# interpolation version
"""
    local_binary_pattern([f], X, npoints, radius, interpolation=Linear(); [rotation], [uniform_degree])

Compute the local binary pattern of gray image `X` using the interpolation-based original
method with circular neighborhood matrix.

This produces better result for `rotation=true` case but is usually slower than the plain
3x3 matrix version `local_binary_pattern(X)`.

# Arguments

- `npoints::Int`(4 ≤ npoints ≤ 32): the number of (uniform-spaced) neighborhood points.
- `radius::Real`(radius ≥ 1.0): the radius of the circular. Larger radius computes the
    pattern of a larger local window/block.
- `interpolation::Union{Degree, InterpolationType}=Linear()`: the interpolation method used
    to generate non-grid pixel value. In most cases, `Linear()` are good enough. One can
    also try other costly interpolation methods, e.g., `Cubic(Line(OnGrid()))`(also known as
    "bicubic"), `Lanczos()`. See also Interpolations.jl for more choices.

!!! info "parameter choices"
    The following parameters are used in [1], with `interpolation=Linear()`.

    | `npoints` | `radius` |
    | --- | --- |
    | ``4`` | ``1.0`` |
    | ``8`` | ``1.0`` |
    | ``12`` | ``1.5`` |
    | ``16`` | ``2.0`` |
    | ``24`` | ``3.0`` |

!!! note "neighborhood order differences"
    Different implementation might use different neighborhood orders; this will change the
    encoding result but will not change the overall distribution. For instance,
    `local_binary_pattern(X)` differs from `local_binary_pattern(X, 8, 1, Constant())` only
    by how `offsets` (see below) are ordered; the former uses column-major top-left to
    bottom-right 3x3 matrix order and the latter uses circular order.

# Examples

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> X = [6 7 9; 5 6 3; 2 1 7]
3×3 $(Matrix{Int}):
 6  7  9
 5  6  3
 2  1  7

julia> local_binary_pattern(X, 4, 1) # 4-neighbor with circular radius 1
3×3 $(Matrix{UInt32}):
 0x00000001  0x00000001  0x00000000
 0x00000003  0x00000002  0x0000000e
 0x00000002  0x00000007  0x00000000

julia> local_binary_pattern(X, 4, 1; rotation=true)
3×3 $(Matrix{UInt32}):
 0x00000001  0x00000001  0x00000000
 0x00000003  0x00000001  0x00000007
 0x00000001  0x00000007  0x00000000
```

# References

- [1] T. Ojala, M. Pietikainen, and T. Maenpaa, “Multiresolution gray-scale and rotation invariant texture classification with local binary patterns,” _IEEE Trans. Pattern Anal. Machine Intell._, vol. 24, no. 7, pp. 971–987, Jul. 2002, doi: 10.1109/TPAMI.2002.1017623.
"""
function local_binary_pattern(
        f, X::AbstractArray, npoints::Int, radius::Real, interpolation=Linear();
        rotation::Bool=false, uniform_degree::Union{Nothing,Int}=nothing)
    interpolation = wrap_BSpline(interpolation)
    offsets = _circular_neighbor_offsets(npoints, radius)
    if interpolation == BSpline(Constant())
        # performance optimization: skip interpolation when it is nearest neighbor interpolation
        offsets = map(offsets) do o
            round.(Int, o, RoundToZero)
        end
        itp = X
    else
        itp = interpolate(X, interpolation)
    end
    # For the sake of type-stability, hardcode to UInt32 here.
    lookups = _build_local_binary_pattern_lookup(UInt32, npoints, rotation=rotation, uniform_degree=uniform_degree)
    local_binary_pattern!(similar(X, UInt32), LBP{ndims(X)}(f), itp, offsets, lookups)
end
function local_binary_pattern(X::AbstractArray, npoints::Int, radius::Real, interpolation=Linear(); kwargs...)
    local_binary_pattern(_default_lbp_mode, X, npoints, radius, interpolation; kwargs...)
end

# Multi-blocks version
"""
    local_binary_pattern([f], X, block_size; [rotation], [uniform_degree])

Compute the local binary pattern of gray image `X` of image blocks with size `block_size`.
This method is usually called MB-LBP[1] or Locally Assembled Binary (LAB) feature[2] in the
literature.

# Arguments

- `X::AbstractMatrix`: the input image matrix. For colorful images, one can manually convert
  it to some monochrome space, e.g., `Gray`, the L-channel of `Lab`. One could also do
  channelwise LBP and then concatenate together.
- `block_size::Tuple{Int,Int}`: positive integers that used to specify the size of each
  block. When `block_size == (1, 1)`, this method degenerates to the pixel-version
  [`local_binary_pattern`](@ref).

For the meaning and values of parameters `rotation` and `uniform_degree` please see the docs
of [`local_binary_pattern`](@ref).

# Examples

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> X = zeros(Int, 9, 9); X[1:3, 1:3] .= 50; X[4:6, 4:6] .= 1; X[7:9, 7:9] .= 50; X
9×9 $(Matrix{Int}):
 50  50  50  0  0  0   0   0   0
 50  50  50  0  0  0   0   0   0
 50  50  50  0  0  0   0   0   0
  0   0   0  1  1  1   0   0   0
  0   0   0  1  1  1   0   0   0
  0   0   0  1  1  1   0   0   0
  0   0   0  0  0  0  50  50  50
  0   0   0  0  0  0  50  50  50
  0   0   0  0  0  0  50  50  50

julia> local_binary_pattern(X, (3, 3))[4, 4] # 0b1000_0001
0x81

julia> local_binary_pattern(X, (3, 3); rotation=true)[4, 4] # 0b0000_0011
0x03

julia> local_binary_pattern(X, (3, 3); uniform_degree=2)[4, 4] # 0x09 (i.e., miscellaneous pattern)
0x09
```

# Extended help

The following is how the block-version of local binary pattern computed with `block_size=(3,
3)`. Check out [`local_binary_pattern`](@ref) for more details of the original pixel-version of
local binary pattern.

```text
2  3  2  |  7  6  6  |  2  1  3
1  2  3  |  5  3  4  |  7  7  3
6  5  4  |  1  4  5  |  1  5  8
-------------------------------     sum                compare          weighted          sum
1  2  3  |  1  7  3  |  8  3  7            28  41  37           0  1  0           0 8 0
1  2  1  |  4  3  6  |  5  2  2     ==>    20  38  45    ==>    0  x  1     ==>   0 x 64  ==> 200
1  2  7  |  5  3  6  |  4  7  7            33  31  39           0  0  1           0 0 128
-------------------------------
4  2  4  |  8  1  5  |  3  7  7
5  7  4  |  1  1  1  |  3  4  3
3  3  1  |  5  1  8  |  1  8  3
```

# References

- [1] Zhang, Lun, et al. "Face detection based on multi-block lbp representation." _International conference on biometrics_. Springer, Berlin, Heidelberg, 2007.
- [2] Yan, Shengye, et al. "Locally assembled binary (LAB) feature with feature-centric cascade for fast and accurate face detection." _2008 IEEE Conference on Computer Vision and Pattern Recognition_. IEEE, 2008.
"""
function local_binary_pattern(f, X::AbstractMatrix, block_size::Dims{2}; rotation::Bool=false, uniform_degree::Union{Nothing,Int}=nothing)
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    lookups = _build_local_binary_pattern_lookup(UInt8, 8; rotation=rotation, uniform_degree=uniform_degree)
    local_binary_pattern!(similar(X, UInt8), LBP(f, block_size), X, offsets, lookups)
end
local_binary_pattern(X::AbstractMatrix, block_size::Dims{2}; kwargs...) = local_binary_pattern(_default_lbp_mode, X, block_size; kwargs...)

function local_binary_pattern!(out::AbstractMatrix, lbp::LocalBinaryPattern, X::AbstractMatrix, offsets, lookups)
    length(offsets) > 32 && throw(ArgumentError("length(offsets)=$(length(offsets)) is expected to be no larger than 32."))

    block_size = lbp.block_size
    offsets = map(offsets) do o
        o .* block_size
    end
    bo = CartesianIndex(block_size .- 1)

    outerR = CartesianIndices(X)
    r = CartesianIndex(
        ceil(Int, maximum(abs.(extrema(first.(offsets))))),
        ceil(Int, maximum(abs.(extrema(last.(offsets)))))
    )
    innerR = first(outerR)+r:last(outerR)-r-bo

    if block_size == (1, 1)
        # TODO(johnnychen94): use LoopVectorization
        @inbounds @simd for I in innerR
            out[I] = lbp_encode(lbp, X, I, offsets)
        end
        for I in EdgeIterator(outerR, innerR)
            out[I] = lbp_encode(lbp, X, I, offsets)
        end
    else
        # It's computationally more efficient to use integral array to compute
        # the mean block value.
        iX = IntegralArray(X)
        # TODO(johnnychen94): it can be faster if we can enable SIMD here.
        #                     Need to skip the boundary check in IntegralArray.
        @inbounds for I in innerR
            out[I] = mlbp_encode_inbounds(lbp, iX, I, bo, offsets)
        end
        for I in EdgeIterator(outerR, innerR)
            out[I] = mlbp_encode(lbp, iX, I, bo, offsets)
        end
    end

    for F in lookups
        @. out = F[out]
    end
    return out
end


###
# Utils
###

function _build_local_binary_pattern_lookup(
        ::Type{T}, nbits;
        rotation::Bool,
        uniform_degree::Union{Nothing,Int}
    ) where T <: Unsigned
    # TODO(johnnychen94): fuse these operation into one broadcasting could potentially give
    # better performance.
    lookups = []
    if rotation
        push!(lookups, rotation_encoding_table(T, nbits))
    end
    if !isnothing(uniform_degree)
        push!(lookups, uniform_encoding_table(T, nbits, uniform_degree))
    end
    return lookups
end

"""
    average_mode(X, I, offsets) -> value
    local_binary_pattern(average_mode, X, args...; kwargs...)

Compute the local binary pattern using the average mode of the block.

Original local binary pattern compares the neighbors with the center value `X[I]`, this
modified version instead uses the mean value of the block.

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> X = [6 7 9; 5 6 3; 2 1 7]
3×3 $(Matrix{Int}):
 6  7  9
 5  6  3
 2  1  7

julia> local_binary_pattern(average_mode, X)
3×3 $(Matrix{UInt8}):
 0xc0  0x40  0x00
 0x68  0xa9  0x1b
 0x28  0x69  0x01
```

See also [`local_binary_pattern`](@ref) for additional arguments and parameters.
"""
Base.@propagate_inbounds function average_mode(X, I::CartesianIndex, offsets)
    v = X[I]
    ax = axes(X)
    rst = v + mapreduce(+, offsets) do o
        p = I.I .+ o
        inbounds = map(in, p, ax)
        all(inbounds) ? _inbounds_getindex(X, p) : v
    end
    return rst / (length(offsets) + 1)
end

Base.@propagate_inbounds function average_mode(X::IntegralArray, I::CartesianIndex, bo::CartesianIndex, offsets)
    # multi-block version
    R = CartesianIndices(X)
    Rfirst, Rlast = first(R), last(R)
    Rp = I:min(I+bo, Rlast)
    v = X[I..last(Rp)]/length(Rp)
    rst = v + mapreduce(+, offsets) do o
        p = CartesianIndex(I.I .+ o)
        Rp = max(p, Rfirst):min(p+bo, Rlast)
        !isempty(Rp) ? X[first(Rp)..last(Rp)]/length(Rp) : v
    end
    return rst / (length(offsets) + 1)
end
