"""
    center_mode(X::AbstractArray, I::CartesianIndex, offsets)
    center_mode(X::IntegralArray, I::CartesianIndex, bo::CartesianIndex, offsets)

Compute the local binary pattern using center pixel as the block mode. This is the
version that original LBP is used.

See also [`average_mode`](@ref) for mLBP (modified LBP) version.
```
"""
Base.@propagate_inbounds center_mode(X, I, offsets) = X[I]
Base.@propagate_inbounds function center_mode(X::IntegralArray, I, bo, offsets)
    Rp = I:min(I+bo, last(CartesianIndices(X)))
    X[first(Rp)..last(Rp)]/length(Rp)
end


"""
    average_mode(X, I, offsets) -> value
    local_binary_pattern(average_mode, X, args...; kwargs...)

Compute the local binary pattern using the average mode of the block. This is also called
mLBP (modified LBP) in the literature [1].

Original local binary pattern compares the neighbors with the center value `X[I]`, this
modified version instead uses the mean value of the block.

See also [`center_mode`](@ref) for original LBP version.

# References

- [1] Wikipedia contributors. "Local binary patterns." _Wikipedia, The Free Encyclopedia_. Wikipedia, The Free Encyclopedia, 17 Nov. 2021. Web. 1 Jan. 2022.
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
