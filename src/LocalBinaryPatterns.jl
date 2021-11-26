module LocalBinaryPatterns

using TiledIteration: EdgeIterator
using ImageCore: GenericGrayImage, OffsetArray

export lbp_original, lbp_original!

include("lbp_original.jl")

end
