using LocalBinaryPatterns: uniform_encoding_table

@testset "uniform_encoding" begin
    lookup = uniform_encoding_table(UInt8, 0:255, 2)
    new_lookup = uniform_encoding_table(UInt8, 0:255, 2)

    # repeat calls are cached in global dictionary
    @test lookup === new_lookup

    X_UInt8 = uniform_encoding_table(UInt8, 0:255, 2)
    @test X_UInt8[0b0000_0001] == 0b0000_0001
    @test X_UInt8[0b0000_0101] == 9

    X_UInt8 = uniform_encoding_table(UInt8, 0:255, 4)
    @test X_UInt8[0b0000_0001] == 0b0000_0001
    @test X_UInt8[0b0000_0101] == 0b0000_0101
    @test X_UInt8[0b0000_1101] == 0b0000_1101
    @test X_UInt8[0b0100_1101] == 9

    X_UInt16 = uniform_encoding_table(UInt16, 0:255, 2)
    @test X_UInt16[0b0000_0001] == 0b0000_0001
    @test X_UInt16[0b0000_0101] == 17
end
