@testset "lbp_rotation_invariant" begin
    # reference result comes from [1]
    # - [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [3 1 0; 13 53 27; 5 91 0]

    out = lbp_rotation_invariant(X)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    fill!(out, 0)
    lbp_rotation_invariant!(out, X)
    @test out == ref_out

    # intensity-invariant
    @test lbp_rotation_invariant(X./255) == out
    # Gray inputs
    @test lbp_rotation_invariant(Gray.(X./255)) == out

    # LBP is only defined for scalar values
    @test_throws MethodError lbp_rotation_invariant(RGB.(Gray.(X./255)))

    # not yet ready for N-dimensional array (although it's doable)
    @test_throws MethodError lbp_rotation_invariant(rand(3, 3, 3))

    @testset "OffsetArrays" begin
        Xo = OffsetArray(X, -1, -1)
        out = lbp_rotation_invariant(Xo)
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == ref_out
    end
end
