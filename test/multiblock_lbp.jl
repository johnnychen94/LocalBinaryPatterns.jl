@testset "multiblock_lbp" begin
    # Only the center-pixel [4, 4] is well-defined in the reference paper [1].
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
        0    0    0  214  214  214  22  22  22
        0    0  128  146  150  150  22  22  22
        0  128  128  130  146  150  22  22  22
      248  200  136  129  145  147  31  31  31
      248  232  200  193  209  209  16  20  20
      248  232  232  201  209  208  16  16  20
      104  104  104  107   64   64   0   0   0
      104  104  104  107   96   64   0   0   0
      104  104  104  107   96   96   0   0   0
    ]

    out = multiblock_lbp(X, (3, 3))
    @test eltype(out) == UInt8
    @test size(out) == (9, 9)
    @test out == ref_out

    ref_out = [
        0    0    6  214  214  214  214  22  22
        0    0    2  214  214  214  214  22  22
       40    8   11  147   23   22  255  31  31
      248  248  201  129  129  150  255  31  31
      248  248  105  129  129  147  255  31  31
      248  248  104  232  201  208  208  20  22
      248  248  255  255  255  208  208  16  22
      104  104  107  107  107   96   64   0   2
      104  104  107  107  107  104  104   8  11
    ]
    @test ref_out == multiblock_lbp(X, (2, 2))

    # ensure default parameter values are not changed accidently
    @test out == multiblock_lbp(X, (3, 3); rotation=false, uniform_degree=nothing)

    # intensity-invariant
    @test multiblock_lbp(X./255, (3, 3)) == out
    # Gray inputs
    @test multiblock_lbp(Gray.(X./255), (3, 3)) == out

    # when `block_size == (1, 1)`, it degenerates to the original pixel version
    @test multiblock_lbp(X, (1, 1)) == local_binary_pattern(X)

    # LBP is only defined for scalar values
    @test_throws MethodError multiblock_lbp(RGB.(Gray.(X./255)), (3, 3))

    # check block_size
    @test_throws ArgumentError multiblock_lbp(X, (0, 3))
    @test_throws MethodError multiblock_lbp(X, (1.5, 1.5))
    @test_throws MethodError multiblock_lbp(X, 3)
    @test_throws MethodError multiblock_lbp(X, (3, ))

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
            0   0   0  91  91  91  11  11  11
            0   0   1  37  45  45  11  11  11
            0   1   1   5  37  45  11  11  11
           31  25  17   3  25  39  31  31  31
           31  29  25   7  29  29   1   5   5
           31  29  29  39  29  13   1   1   5
           13  13  13  91   1   1   0   0   0
           13  13  13  91   3   1   0   0   0
           13  13  13  91   3   3   0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=true)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end

    @testset "Uniform encoding" begin
        ref_out = [
            0    0    0  9   9   9   9   9   9
            0    0  128  9   9   9   9   9   9
            0  128  128  9   9   9   9   9   9
          248    9    9  9   9   9  31  31  31
          248    9    9  9   9   9  16   9   9
          248    9    9  9   9   9  16  16   9
            9    9    9  9  64  64   0   0   0
            9    9    9  9  96  64   0   0   0
            9    9    9  9  96  96   0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=false, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end

    @testset "Rotation Invariant, Uniform encoding" begin
        ref_out = [
            0  0  0  9  9  9   9   9   9
            0  0  1  9  9  9   9   9   9
            0  1  1  9  9  9   9   9   9
           31  9  9  3  9  9  31  31  31
           31  9  9  7  9  9   1   9   9
           31  9  9  9  9  9   1   1   9
            9  9  9  9  1  1   0   0   0
            9  9  9  9  3  1   0   0   0
            9  9  9  9  3  3   0   0   0
        ]

        out = multiblock_lbp(X, (3, 3); rotation=true, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end
end
