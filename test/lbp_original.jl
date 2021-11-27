@testset "lbp_original" begin
    # reference result comes from [1]
    # - [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [192 64 0; 104 169 27; 40 107 0]

    out = lbp_original(X)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    # ensure default parameter values are not changed accidently
    @test out == lbp_original(X; rotation=false)

    fill!(out, 0)
    lbp_original!(out, X)
    @test out == ref_out

    # intensity-invariant
    @test lbp_original(X./255) == out
    # Gray inputs
    @test lbp_original(Gray.(X./255)) == out

    # LBP is only defined for scalar values
    @test_throws MethodError lbp_original(RGB.(Gray.(X./255)))

    # not yet ready for N-dimensional array (although it's doable)
    @test_throws MethodError lbp_original(rand(3, 3, 3))

    @testset "OffsetArrays" begin
        Xo = OffsetArray(X, -1, -1)
        out = lbp_original(Xo)
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == ref_out
    end

    @testset "Rotation Invariant" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [3 1 0; 13 53 27; 5 91 0]

        out = lbp_original(X; rotation=true)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [192 64 0; 9 9 9; 9 9 0]

        out = lbp_original(X; rotation=false, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Rotation Invariant, Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [3 1 0; 9 9 9; 9 9 0]

        out = lbp_original(X; rotation=true, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end
end
