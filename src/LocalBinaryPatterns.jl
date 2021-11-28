module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using Interpolations
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays

export lbp_original

include("lbp_original.jl")

include("bitrotate_encoding.jl")
include("uniform_encoding.jl")
include("utils.jl")
include("compat.jl")

end
