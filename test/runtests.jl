using LocalBinaryPatterns
using ImageCore
using ImageCore.OffsetArrays
using IntegralArrays
using Interpolations
using Statistics
using Test
using Aqua, Documenter

@testset "LocalBinaryPatterns.jl" begin

@testset "meta-quality" begin
    Aqua.test_ambiguities(LocalBinaryPatterns)
    Aqua.test_all(LocalBinaryPatterns; ambiguities=false)
    if VERSION <= v"1.8.0-DEV.1070"
        # TODO(johnnychen94): remove this version check when `@_inline_meta` deprecation is handled in upstream ecosystems
        doctest(LocalBinaryPatterns; manual = false)
    end
end

include("bitrotate_encoding.jl")
include("uniform_encoding.jl")
include("modes.jl")
include("local_binary_pattern.jl")

end
