using LocalBinaryPatterns
using ImageCore
using Test
using Aqua, Documenter

@testset "LocalBinaryPatterns.jl" begin

@testset "meta-quality" begin
    Aqua.test_ambiguities(LocalBinaryPatterns)
    Aqua.test_all(LocalBinaryPatterns; ambiguities=false)
    doctest(LocalBinaryPatterns; manual = false)
end

include("lbp_original.jl")

end
