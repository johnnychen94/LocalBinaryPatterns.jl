using LocalBinaryPatterns
using Documenter

format = Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    assets=String[],
)
makedocs(;
    sitename="LocalBinaryPatterns Documentation",
    modules=[LocalBinaryPatterns],
    format=format,
    pages=[
        "Home" => "index.md",
        "API References" => "reference.md"
    ],
)

deploydocs(;repo="github.com/johnnychen94/LocalBinaryPatterns.jl")
