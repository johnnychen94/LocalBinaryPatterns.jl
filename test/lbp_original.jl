@testset "lbp_original" begin
    # reference result comes from [1]
    # - [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [192 64 0; 104 169 27; 40 107 0]

    out = lbp_original(X)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    fill!(out, 0)
    lbp_original!(out, X)
    @test out == ref_out

    # intensity-invariant
    @test lbp_original(X./255) == out
    # Gray inputs
    @test lbp_original(Gray.(X./255)) == out

    # LBP is only defined for scalar values
    @test_throws MethodError lbp_original(RGB.(Gray.(X./255)))
end
