"""
    multiblock_lbp(X, block_size; [rotation], [uniform_degree])

Compute the local binary pattern of gray image `X` of image blocks with size `block_size`.
This method is usually called MB-LBP[1] or Locally Assembled Binary (LAB) feature[2] in the
literature.

# Arguments

- `X::AbstractMatrix`: the input image matrix. For colorful images, one can manually convert
  it to some monochrome space, e.g., `Gray`, the L-channel of `Lab`. One could also do
  channelwise LBP and then concatenate together.
- `block_size::Tuple{Int,Int}`: positive integers that used to specify the size of each
  block. When `block_size == (1, 1)`, this method degenerates to the pixel-version
  [`lbp_original`](@ref).

For the meaning and values of parameters `rotation` and `uniform_degree` please see the docs
of [`lbp_original`](@ref).

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

julia> multiblock_lbp(X, (3, 3))[4, 4] # 0b1000_0001
0x81

julia> multiblock_lbp(X, (3, 3); rotation=true)[4, 4] # 0b0000_0011
0x03

julia> multiblock_lbp(X, (3, 3); uniform_degree=2)[4, 4] # 0x09 (i.e., miscellaneous pattern)
0x09
```

# Extended help

The following is how the block-version of local binary pattern computed with `block_size=(3,
3)`. Check out [`lbp_original`](@ref) for more details of the original pixel-version of
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
function multiblock_lbp(X::AbstractArray, block_size::Dims{2}; rotation::Bool=false, uniform_degree::Union{Nothing,Int}=nothing)
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    lookups = _build_lbp_original_lookup(UInt8, 8; rotation=rotation, uniform_degree=uniform_degree)
    multiblock_lbp!(similar(X, UInt8), X, offsets, lookups, block_size)
end

function multiblock_lbp!(out, X::AbstractMatrix, offsets, lookups, block_size::Dims{2})
    all(block_size .> 0) || throw(ArgumentError("block size `block_size=$(block_size)` should be positive numbers."))
    offsets = map(offsets) do o
        o .* block_size
    end

    # It's computationally more efficient to use integral array to compute
    # the mean block value.
    iX = IntegralArray(X)

    outerR = CartesianIndices(X)
    r = CartesianIndex(
        ceil(Int, maximum(abs.(extrema(first.(offsets))))),
        ceil(Int, maximum(abs.(extrema(last.(offsets)))))
    )

    bo = CartesianIndex(block_size .- 1)
    innerR = first(outerR)+r:last(outerR)-r-bo

    # TODO(johnnychen94): use LoopVectorization
    @inbounds for I in innerR
        # since the block is regular for inner region, doing sum is equivalent to doing average.
        gc = iX[I..I+bo]
        # This inner loop fuses the binary pattern build stage and encoding stage for
        # better performance.
        rst = 0
        # TODO(johnnychen94): it can be faster if we can enable SIMD here.
        #                     Need to skip the boundary check in IntegralArray.
        for i in 1:length(offsets)
            p = CartesianIndex(I.I .+ offsets[i])
            rst += ifelse(gc <= iX[p..p+bo], 1, 0) << (i-1)
        end
        out[I] = rst
    end

    # boundary conditions are undefined in [1]; here we directly skip out-of-boundary values
    @inbounds for I in EdgeIterator(outerR, innerR)
        R = I:min(I+bo, last(outerR))
        gc = iX[first(R)..last(R)]/length(R)
        rst = 0
        for i in 1:length(offsets)
            p = CartesianIndex(I.I .+ offsets[i])
            checkbounds(Bool, X, p.I...) || continue
            Rp = p:min(p+bo, last(outerR))
            gp = iX[first(Rp)..last(Rp)]/length(Rp)
            rst += ifelse(gc <=gp, 1, 0) << (i-1)
        end
        out[I] = rst
    end

    for F in lookups
        @. out = F[out]
    end

    return out
end
