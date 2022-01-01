@testset "modes" begin

@testset "center_mode" begin
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    X = [6  7  9
         5  6  2
         2  1  6]
    @test 6 == center_mode(X, CartesianIndex(2, 2), offsets)

    iX = IntegralArray(X)
    @test 6 == center_mode(iX, CartesianIndex(2, 2), CartesianIndex(0, 0), offsets)
    @test mean(X[2:3,2:3]) == center_mode(iX, CartesianIndex(2, 2), CartesianIndex(1, 1), offsets)
end

@testset "average_mode" begin
    offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))
    X = [6  7  9
         5  6  2
         2  1  6]
    @test mean(X) ≈ average_mode(X, CartesianIndex(2, 2), offsets)

    iX = IntegralArray(X)
    @test mean(X) ≈ average_mode(iX, CartesianIndex(2, 2), CartesianIndex(0, 0), offsets)

    ref = (mean(X[1:2, 1:2]) + mean(X[2:3, 1:2]) + mean(X[3, 1:2]) +
          mean(X[1:2, 2:3]) + mean(X[2:3, 2:3]) + mean(X[3, 2:3]) +
          mean(X[1:2, 3  ]) + mean(X[2:3, 3  ]) + mean(X[3, 3]))/9
    @test ref ≈ average_mode(iX, CartesianIndex(2, 2), CartesianIndex(1, 1), offsets)
end

@testset "LBP" begin
    lbp = LBP((1, 1))
    @test lbp.mode == center_mode
    @test lbp.block_size == (1, 1)

    @test_throws MethodError LBP()
    @test_throws ArgumentError LBP((0, 3))
end

end
