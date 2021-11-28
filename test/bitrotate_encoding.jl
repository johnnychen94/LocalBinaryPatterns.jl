using LocalBinaryPatterns: RotationEncodingTable, rotation_encoding_table

@testset "bitrotate_encoding" begin
    @testset "RotationEncodingTable" begin
        X_UInt8 = RotationEncodingTable{UInt8}(0:255)
        X_UInt16 = RotationEncodingTable{UInt16}(0:255)

        @test length(unique(X_UInt8)) == 36
        @test length(unique(X_UInt16)) == 129
        @test X_UInt8[0b1010_1001] == 0b0011_0101
        @test X_UInt8[0b1111_1101] != X_UInt16[0b0000_0000_1111_1101]

        X_UInt8_half = RotationEncodingTable{UInt8}(128:255)
        @test all(128 .<= X_UInt8_half .<= 255)
        @test length(unique(X_UInt8_half)) == 35
        @test X_UInt8_half[0b1010_1001] == 0b1001_1010

        @test_throws ArgumentError RotationEncodingTable{UInt8}(0:65535)

        X = RotationEncodingTable{UInt32}(0:2^24-1) # 64MB storage for 24 bits
        # ensure the getindex calculation finishes in finite time
        t = @elapsed sum(X)
        @test t < 10 # is expected to be smaller than 1

        # everytime it allocates new memory to store the result
        @test RotationEncodingTable{UInt8}(0:255) !== RotationEncodingTable{UInt8}(0:255)
    end

    # repeat calls are cached in global dictionary
    lookup = rotation_encoding_table(UInt8, 0:255)
    new_lookup = rotation_encoding_table(UInt8, 0:255)
    @test lookup === new_lookup
    @test lookup === rotation_encoding_table(UInt8, 8)
end
