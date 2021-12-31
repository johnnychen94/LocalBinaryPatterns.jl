module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using Interpolations
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays
using IntegralArrays

export local_binary_pattern

include("local_binary_pattern.jl")

include("bitrotate_encoding.jl")
include("uniform_encoding.jl")
include("utils.jl")
include("compat.jl")

end
