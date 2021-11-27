module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays

export
    lbp_original, lbp_original!,
    lbp_ri, lbp_ri!

include("lbp_original.jl")
include("lbp_ri.jl")
include("utils.jl")

end
