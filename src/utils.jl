function _circular_neighbor_offsets(npoints::Integer, radius::Real)
    radius >= 1 || throw(ArgumentError("radius >= 1.0 is expected."))

    ntuple(npoints) do p
        θ=2π*(p-1)/npoints
        pos = (-1, 1) .* radius .* sincos(θ)
        # rounding the positions for more stable result when applying interpolation
        round.(pos; digits=8)
    end
end

# A helper function to let `interpolation` keyword correctly understands `Degree` inputs
@inline wrap_BSpline(itp::Interpolations.InterpolationType) = itp
@inline wrap_BSpline(degree::Interpolations.Degree) = BSpline(degree)

@inline _inbounds_getindex(A, I) = @inbounds A[I...]
@inline _inbounds_getindex(A::Interpolations.AbstractInterpolation, I) = @inbounds A(I...)
