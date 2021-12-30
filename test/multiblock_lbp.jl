@testset "multiblock_lbp" begin
    # Only the center-value is well-defined in the reference paper [1].
    # [1] Zhang, Lun, et al. "Face detection based on multi-block lbp representation." _International conference on biometrics_. Springer, Berlin, Heidelberg, 2007.
    #
    # 50  0  0      1 0 0    1 0 0
    #  0  1  0   => 0 x 0 => 0 x 0   => 0b10000001 (129)
    #  0  0  50     0 0 1    0 0 128
    X = [
        50  50  50  0  0  0   0   0   0
        50  50  50  0  0  0   0   0   0
        50  50  50  0  0  0   0   0   0
         0   0   0  1  1  1   0   0   0
         0   0   0  1  1  1   0   0   0
         0   0   0  1  1  1   0   0   0
         0   0   0  0  0  0  50  50  50
         0   0   0  0  0  0  50  50  50
         0   0   0  0  0  0  50  50  50
    ]
    ref_out = [
        0    0    0    6  214  214  22  22  22
        0    0    0    2  214  214  22  22  22
        0    0    0  130  146  150  22  22  22
       40    8  136  139  131  147  23  31  31
      248  248  200  137  129  145  19  31  31
      248  248  232  201  193  209  17  16  20
      104  104  104  105   73   65   0   0   0
      104  104  104  107  107   64   0   0   0
      104  104  104  107  107   96   0   0   0
    ]

    out = multiblock_lbp(X, (3, 3))
    @test eltype(out) == UInt8
    @test size(out) == (9, 9)
    @test out == ref_out

    # ensure default parameter values are not changed accidently
    @test out == multiblock_lbp(X, (3, 3); rotation=false, uniform_degree=nothing)

    # intensity-invariant
    @test multiblock_lbp(X./255, (3, 3)) == out
    # Gray inputs
    @test multiblock_lbp(Gray.(X./255), (3, 3)) == out

    # when `block_size == (1, 1)`, it degenerates to the original pixel version
    @test multiblock_lbp(X, (1, 1)) == lbp_original(X)

    # LBP is only defined for scalar values
    @test_throws MethodError multiblock_lbp(RGB.(Gray.(X./255)), (3, 3))

    # not yet ready for N-dimensional array (although it's doable)
    @test_throws MethodError multiblock_lbp(rand(3, 3, 3), (3, 3))

    @testset "OffsetArrays" begin
        Xo = OffsetArray(X, -1, -1)
        out = multiblock_lbp(Xo, (3, 3))
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == ref_out
    end

    @testset "Rotation Invariant" begin
        ref_out = [
            0   0   0   3  91  91  11  11  11
            0   0   0   1  91  91  11  11  11
            0   0   0   5  37  45  11  11  11
            5   1  17  23   7  39  23  31  31
           31  31  25  19   3  25  19  31  31
           31  31  29  39   7  29  17   1   5
           13  13  13  45  37   5   0   0   0
           13  13  13  91  91   1   0   0   0
           13  13  13  91  91   3   0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=true)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end

    @testset "Uniform encoding" begin
        ref_out = [
            0    0  0  6  9   9  9   9   9
            0    0  0  2  9   9  9   9   9
            0    0  0  9  9   9  9   9   9
            9    8  9  9  9   9  9  31  31
          248  248  9  9  9   9  9  31  31
          248  248  9  9  9   9  9  16   9
            9    9  9  9  9   9  0   0   0
            9    9  9  9  9  64  0   0   0
            9    9  9  9  9  96  0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=false, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end

    @testset "Rotation Invariant, Uniform encoding" begin
        ref_out = [
            0   0  0  3  9  9  9   9   9
            0   0  0  1  9  9  9   9   9
            0   0  0  9  9  9  9   9   9
            9   1  9  9  7  9  9  31  31
           31  31  9  9  3  9  9  31  31
           31  31  9  9  7  9  9   1   9
            9   9  9  9  9  9  0   0   0
            9   9  9  9  9  1  0   0   0
            9   9  9  9  9  3  0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=true, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end
end
