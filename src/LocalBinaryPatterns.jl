module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using Interpolations
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays
using IntegralArrays

export
    local_binary_pattern,
    LBP,
    center_mode,
    average_mode

include("local_binary_pattern.jl")

include("block_modes.jl")
include("bitrotate_encoding.jl")
include("uniform_encoding.jl")
include("utils.jl")
include("compat.jl")

end
