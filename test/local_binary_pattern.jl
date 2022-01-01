using LocalBinaryPatterns: average_mode

@testset "local_binary_pattern" begin
    # reference result comes from [1]
    # - [1] T. Ojala, M. Pietikäinen, and D. Harwood, “A comparative study of texture measures with classification based on featured distributions,” _Pattern Recognition_, vol. 29, no. 1, pp. 51–59, Jan. 1996, doi: 10.1016/0031-3203(95)00067-4.
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [192 64 0; 104 169 27; 40 107 0]

    out = local_binary_pattern(X)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    # ensure default parameter values are not changed accidently
    @test out == local_binary_pattern(X; rotation=false, uniform_degree=nothing)

    # intensity-invariant
    @test local_binary_pattern(X./255) == out
    # Gray inputs
    @test local_binary_pattern(Gray.(X./255)) == out

    # LBP is only defined for scalar values
    @test_throws MethodError local_binary_pattern(RGB.(Gray.(X./255)))

    # not yet ready for N-dimensional array (although it's doable)
    @test_throws MethodError local_binary_pattern(rand(3, 3, 3))

    @testset "OffsetArrays" begin
        Xo = OffsetArray(X, -1, -1)
        out = local_binary_pattern(Xo)
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == ref_out
    end

    @testset "Rotation Invariant" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [3 1 0; 13 53 27; 5 91 0]

        out = local_binary_pattern(X; rotation=true)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [192 64 0; 9 9 9; 9 9 0]

        out = local_binary_pattern(X; rotation=false, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Rotation Invariant, Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [3 1 0; 9 9 9; 9 9 0]

        out = local_binary_pattern(X; rotation=true, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (3, 3)
        @test out == ref_out
    end
end

@testset "local_binary_pattern, Interpolation-based" begin
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [1 1 0; 3 2 14; 2 7 0]
    out = local_binary_pattern(X, 4, 1)
    @test eltype(out) == UInt32
    @test size(out) == (3, 3)
    @test out == ref_out

    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [129 1 0; 7 14 124; 6 31 0]
    out = local_binary_pattern(X, 8, 1)
    @test eltype(out) == UInt32
    @test size(out) == (3, 3)
    @test out == ref_out

    @testset "OffsetArrays" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [129 1 0; 7 14 124; 6 31 0]
        Xo = OffsetArray(X, -1, -1)
        out = local_binary_pattern(Xo, 8, 1)
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == ref_out
    end

    # ensure default parameter values are not changed accidently
    out = local_binary_pattern(X, 8, 1)
    @test local_binary_pattern(X, 8, 1, Linear(); rotation=false, uniform_degree=nothing) == out
    @test local_binary_pattern(X, 8, 1, Linear()) == local_binary_pattern(X, 8, 1, BSpline(Linear()))

    # The pattern encoding is different because the offset orders are not the same
    @test local_binary_pattern(X, 8, 1, Constant()) != local_binary_pattern(X)

    @testset "Rotation Invariant" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [129 1 0; 7 7 31; 3 31 0]
        out = local_binary_pattern(X, 8, 1; rotation=true)
        @test eltype(out) == UInt32
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [33 1 0; 33 33 33; 33 33 0]
        out = local_binary_pattern(X, 8, 1; uniform_degree=2)
        @test eltype(out) == UInt32
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "Rotation Invariant, Uniform encoding" begin
        X = [6 7 9; 5 6 3; 2 1 7]
        ref_out = [33 1 0; 33 33 33; 33 33 0]
        out = local_binary_pattern(X, 8, 1; rotation=true, uniform_degree=2)
        @test eltype(out) == UInt32
        @test size(out) == (3, 3)
        @test out == ref_out
    end

    @testset "radius" begin
        X = [3  1  1  2  2
             1  2  1  5  5
             1  1  2  4  5
             3  3  5  1  2
             3  3  5  2  2]
        ref_out = [0   9   13  8   8
                   9   9   13  0   0
                   11  11  9   0   0
                   1   0   0   6   6
                   1   0   0   6   6]
        out = local_binary_pattern(X, 4, 2)
        @test eltype(out) == UInt32
        @test size(out) == (5, 5)
        @test out == ref_out

        @test_nowarn local_binary_pattern(X, 4, 1.5)
        @test_throws ArgumentError local_binary_pattern(X, 4, 0.5)
    end
end

@testset "average mode" begin
    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [192 64 0; 104 169 27; 40 105 1]
    out = local_binary_pattern(LBP(average_mode, (1, 1)), X)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [3 1 0; 13 53 27; 5 45 1]
    out = local_binary_pattern(LBP(average_mode, (1, 1)), X, rotation=true)
    @test eltype(out) == UInt8
    @test size(out) == (3, 3)
    @test out == ref_out

    X = [6 7 9; 5 6 3; 2 1 7]
    ref_out = [1 1 0; 3 6 28; 2 15 0]
    out = local_binary_pattern(LBP(average_mode, (1, 1)), X, 6, 1)
    @test eltype(out) == UInt32
    @test size(out) == (3, 3)
    @test out == ref_out
end

@testset "multiblock LBP" begin
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

    out = local_binary_pattern(LBP((3, 3)), X)
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
    @test ref_out == local_binary_pattern(LBP((2, 2)), X)

    # ensure default parameter values are not changed accidently
    @test out == local_binary_pattern(LBP((3, 3)), X; rotation=false, uniform_degree=nothing)

    # intensity-invariant
    @test_broken local_binary_pattern(LBP((3, 3)), X./3) == out # float-point numerical error
    @test local_binary_pattern(LBP((3, 3)), X./4) == out
    # Gray inputs
    @test local_binary_pattern(LBP((3, 3)), Gray.(X./4)) == out

    # when `block_size == (1, 1)`, it degenerates to the original pixel version
    @test local_binary_pattern(LBP((1, 1)), X) == local_binary_pattern(X)

    # LBP is only defined for scalar values
    @test_throws MethodError local_binary_pattern(LBP((3, 3)), RGB.(Gray.(X./255)))

    # check block_size
    @test_throws ArgumentError local_binary_pattern(LBP((0, 3)), X, (0, 3))
    @test_throws MethodError local_binary_pattern(LBP((1.5, 1.5)), X)
    @test_throws MethodError local_binary_pattern(LBP(3), X)
    @test_throws MethodError local_binary_pattern(LBP((3, )), X)

    # not yet ready for N-dimensional array (although it's doable)
    @test_throws MethodError local_binary_pattern(rand(3, 3, 3), (3, 3))

    @testset "OffsetArrays" begin
        Xo = OffsetArray(X, -1, -1)
        out = local_binary_pattern(LBP((3, 3)), Xo)
        @test axes(out) == axes(Xo)
        @test OffsetArrays.no_offset_view(out) == local_binary_pattern(LBP((3, 3)), X)
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

        out = local_binary_pattern(LBP((3, 3)), X; rotation=true)
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

        out = local_binary_pattern(LBP((3, 3)), X; rotation=false, uniform_degree=2)
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

        out = local_binary_pattern(LBP((3, 3)), X; rotation=true, uniform_degree=2)
        @test eltype(out) == UInt8
        @test size(out) == (9, 9)
        @test out == ref_out
    end

    @testset "Interpolation-based" begin
        @test_throws ArgumentError local_binary_pattern(LBP((3, 3)), X, 4, 1)
        @test_throws ArgumentError local_binary_pattern(LBP((3, 3)), X, 4, 1; rotation=true)
        @test_throws ArgumentError local_binary_pattern(LBP((3, 3)), X, 4, 1; uniform_degree=2)
    end
end

@testset "multiblock LBP + average_mode" begin
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
        0    0    0    2    2    2   4   4   4
        0    0    0  130  130  146  16  20  20
        0    0  128  130  146  144  16  20  20
        8  136  136  129  145  145  16  20  20
        8  136  200  193  209  208  16  20  20
        8  200  192  193  208  208  16  16  20
       32   64   64   64   64   64   0   0   0
       32   96   96   96   96   64   0   0   0
       32   96   96   96   96   96   0   0   0
    ]
    out = local_binary_pattern(LBP(average_mode, (3, 3)), X)
    @test out == ref_out

    # the center pixel is well-defined here
    @test out[4, 4] == local_binary_pattern(LBP(average_mode, (1, 1)), [50 0 0; 0 1 0; 0 0 50])[2, 2]
end
