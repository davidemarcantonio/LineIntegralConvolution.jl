module LineIntegralConvolution

    using Random
    using Plots
    using Images

    include("utils.jl")
    export load_2d_electrostatic_example

    include("functions.jl")
    export lic
    export lic_process

end