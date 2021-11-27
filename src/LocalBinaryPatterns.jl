module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using ImageCore

export lbp_original, lbp_original!

include("lbp_original.jl")

end
