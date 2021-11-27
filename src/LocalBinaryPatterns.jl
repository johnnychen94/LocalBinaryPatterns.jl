module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays

export
    lbp_original, lbp_original!

include("lbp_original.jl")

include("utils.jl")
include("compat.jl")

end
