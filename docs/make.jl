push!(LOAD_PATH, "../src/")

using Documenter
using LineIntegralConvolution

makedocs(
    modules = [LineIntegralConvolution],
    format = Documenter.HTML(),
    sitename="LineIntegralConvolution", 
    authors = "Davide Marcantonio"
)
    