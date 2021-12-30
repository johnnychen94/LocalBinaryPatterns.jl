module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using Interpolations
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays
using IntegralArrays

export lbp_original, multiblock_lbp

include("lbp_original.jl")
include("multiblock_lbp.jl")

include("bitrotate_encoding.jl")
include("uniform_encoding.jl")
include("utils.jl")
include("compat.jl")

end
