"""
    LBP([mode], block_size)

Construct the classic local binary pattern methods. This includes from O-LBP (original LBP)
[1-3], mLBP (modified LBP), MB-LBP (multi-block LBP)[5-6].

# Arguments

- `mode`: a function that defines how block mode value is computed. By default it uses the
center pixel value as the block mode value, a.k.a., the original LBP [1]. Currently
supported modes are: [`center_mode`](@ref), [`average_mode`](@ref).
- `block_size::Dims`: positive integers that used to specify the size of each block [5-6].
When `all(block_size .== 1)`, it degenerates to the pixel-version [1].

# Examples

```jldoctest; setup=:(using LocalBinaryPatterns)
julia> using LocalBinaryPatterns: center_mode, average_mode

julia> X = [6 7 9; 5 6 2; 2 1 6]
3×3 $(Matrix{Int}):
 6  7  9
 5  6  2
 2  1  6

julia> local_binary_pattern(LBP(center_mode, (1, 1)), X)[2, 2] # mode value: 6
0xa9

julia> local_binary_pattern(LBP(average_mode, (1, 1)), X)[2, 2] # mode value: 44//9 ≈ 4.8
0xab
```

# References

- [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
- [2] T. Ojala, M. Pietikäinen, and T. Mäenpää, “A Generalized Local Binary Pattern Operator for Multiresolution Gray Scale and Rotation Invariant Texture Classification,” in _Advances in Pattern Recognition — ICAPR 2001, vol. 2013, S. Singh, N. Murshed, and W. Kropatsch, Eds. Berlin, Heidelberg: Springer Berlin Heidelberg_, 2001, pp. 399–408. doi: 10.1007/3-540-44732-6_41.
- [3] Pietikäinen, Matti, Timo Ojala, and Zelin Xu. "Rotation-invariant texture classification using feature distributions." _Pattern recognition_ 33.1 (2000): 43-52.
- [4] T. Ojala, M. Pietikainen, and T. Maenpaa, “Multiresolution gray-scale and rotation invariant texture classification with local binary patterns,” _IEEE Trans. Pattern Anal. Machine Intell._, vol. 24, no. 7, pp. 971–987, Jul. 2002, doi: 10.1109/TPAMI.2002.1017623.
- [5] Zhang, Lun, et al. "Face detection based on multi-block lbp representation." _International conference on biometrics_. Springer, Berlin, Heidelberg, 2007.
- [6] Yan, Shengye, et al. "Locally assembled binary (LAB) feature with feature-centric cascade for fast and accurate face detection." _2008 IEEE Conference on Computer Vision and Pattern Recognition_. IEEE, 2008.
"""
struct LBP{N,F} <: LocalBinaryPattern
    mode::F # mode(X::AbstractArray{T,N}, I::CartesianIndex{N}, offsets::Tuple)
    block_size::Dims{N}
    function LBP(mode::F, block_size::Dims{N}) where {N,F}
        all(block_size .> 0) || throw(ArgumentError("block size `block_size=$(block_size)` should be positive integers."))
        new{N, F}(mode, block_size)
    end
end
LBP(block_size::Dims) = LBP(center_mode, block_size)

Base.show(io::IO, lbp::LBP) = print(io, "LBP(", lbp.mode, ", ", lbp.block_size, ")")

Base.@propagate_inbounds function _lbp_encode(lbp::LBP, X::AbstractMatrix, I::CartesianIndex{2}, offsets::Tuple)
    rst = 0
    gc = lbp.mode(X, I, offsets)
    for i in 1:length(offsets)
        p = I.I .+ offsets[i]
        @boundscheck checkbounds(Bool, X, p...) || continue
        rst += ifelse(gc <= _inbounds_getindex(X, p), 1, 0) << (i-1)
    end
    return rst
end

Base.@propagate_inbounds function _mlbp_encode_inbounds(
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
function _mlbp_encode(
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
