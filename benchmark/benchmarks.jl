# Usage:
#     julia benchmark/run_benchmarks.jl
using BenchmarkTools
using LocalBinaryPatterns
using ImageTransformations
using TestImages
using ImageCore

include("mapwindow_impl.jl")

on_CI = haskey(ENV, "GITHUB_ACTIONS")

img = testimage("cameraman")
tst_sizes = (256, 512)
tst_types = (N0f8, Float32, Gray{N0f8}, Gray{Float32})

const SUITE = BenchmarkGroup()

alg_list = (( "Original", lbp_original),
            )

function add_algorithm_benchmark!(suite, img, alg_name, alg;
                                  tst_sizes, tst_types)
    haskey(suite, alg_name) || (suite[alg_name] = BenchmarkGroup())

    for T in tst_types
        haskey(suite[alg_name], T) || (suite[alg_name][T] = BenchmarkGroup())
        for sz in tst_sizes
            tst_img = imresize(T.(img), (sz, sz))

            suite[alg_name][T]["$sz×$sz"] = @benchmarkable $alg($tst_img)
        end
    end
end


for (alg_name, alg) in alg_list
    add_algorithm_benchmark!(SUITE, img, alg_name, alg;
                             tst_sizes=tst_sizes,
                             tst_types=tst_types)
end


# compare with the simple mapwindow version
suite = SUITE["Original (mapwindow)"] = BenchmarkGroup()
for T in tst_types
    haskey(suite, T) || (suite[T] = BenchmarkGroup())
    for sz in tst_sizes
        tst_img = imresize(T.(img), (sz, sz))
        suite[T]["$sz×$sz"] = @benchmarkable lbp_origin_mapwindow($tst_img)
    end
end