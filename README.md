# LineIntegralConvolution.jl

Line Integral Convolution (LIC) visualization tool for 2D vector fields

## Documentation

Documentation is available [here](https://davidemarcantonio.github.io/LineIntegralConvolution.jl/)

## Dependencies

This package builds over the following Julia packages

* [Random](https://docs.julialang.org/en/v1/stdlib/Random/) (standard library) 
* [Plots](https://docs.juliaplots.org/latest/tutorial/) v1.38
* [Images](https://juliaimages.org/latest/install/) v0.25

## Example

```julia
using Plots
using LineIntegralConvolution

pt_per_meter = 25  # resolution, number of pixels per meter
x_min = -10.
x_max = 10.
y_min = 2.
y_max = 13.
n_charges = 10  # number of charges to simulate
charge_value = 1e-6  # electric charge value (C)
distribution = "circle"  # "circle", "random"
SEED = 1
field_result = load_2d_electrostatic_example(
    x_min, x_max, y_min, y_max, pt_per_meter,
    n_charges=n_charges, distribution=distribution,
    seed=SEED,
    charge_value=charge_value
)

lic(field_result.f)

savefig("electric_charges.pdf")  # export to PDF vectorial image
savefig("electric_charges.png")  # export to PNG
```

![png](data/electric_charges.png)

## About

This project started on `April, 2023`.
