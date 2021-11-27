module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using ImageCore
using ImageCore.OffsetArrays
using StaticArrays

export
    lbp_original, lbp_original!,
    lbp_rotation_invariant, lbp_rotation_invariant!

include("lbp_original.jl")
include("lbp_rotation_invariant.jl")

include("utils.jl")
include("compat.jl")

end
