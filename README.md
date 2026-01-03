# MultiFlowSolver

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dcunha95.github.io/MultiFlowSolver.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dcunha95.github.io/MultiFlowSolver.jl/dev/)
[![Build Status](https://github.com/dcunha95/MultiFlowSolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/dcunha95/MultiFlowSolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/dcunha95/MultiFlowSolver.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dcunha95/MultiFlowSolver.jl)


MultiFlowSolver is a package designed to solve the Multicommodity Network Flow Problem (MNFP) by a Column Generation Algorithm (CGA)

## Usage

Create a `MnfpData` representing an existing instance, where `A`, `u` and `c` are the instance's arcs, flow capacities and costs: 

```julia
using MultiFlowSolver

vertex_amount::Int64 = 4
A::Vector{Tuple{Int64, Int64}} = [(1,2), (1,3), (2,3), (2,4), (3,4)]
u::Matrix{Float64} = [10 5 15 6 7]
c::Matrix{Float64} = [
    3 2 2 4 1;
    4 2 1 2 5
]
d::Matrix{Float64} = [
    -6 -2 0 8;
    -7 0 3 4
]

mnfp_instance = MnfpData(vertex_amount, A, u, c, d)
```

Call the `mnfp_cga` function to solve it:

```julia
out = mnfp_cga(mnfp_instance, max_time=1200)
out.rmlp_obj
```

## Serialization

The package comes with the functions `save_mnfp`, `save_solution`, `load_mnfp` and `write_mnfp_txt`.

```julia

save_mnfp(mnfp_instance, "instance.jls")
save_solution(out.solution, "instance_solution.jls")

load_mnfp("instance.jls")

write_mnfp_txt("instance.txt", mnfp_instance)
load_mnfp("instance.txt")
```

Functions `save_mnfp`, `save_solution` are simple serialization into `.jls` objects. `write_mnfp_txt` saves the instance exclusively as a text file (note that there is a small data loss in terms of display data). Finally, `load_mnfp` reads both `.jls` files and `.txt` files.

## Random Instance Creation

Alternatively, you may generate a random `MnfpData` instance (Note that generated instances aren't guaranteed to be feasible).

```julia
mnfp_instance = generate_planar_grid_mnfp(
    vertex_amount=30,
    K=20,
    random_capacity_slack=0.5,
    base_capacity_slack=5.0,
)
out = mnfp_cga(mnfp_instance)
out.rmlp_obj
```

If an instance is too easy (or infeasible), `harden_instance` will attempt to increase its difficulty while ensuring feasibility, modifying the instance's capacity and saving it and its picture on `instance_base_path` and `instance_base_img_path` respectively.

```julia
harden_instance(
    mnfp_instance, 
    instance_base_path, 
    instance_base_img_path, 
    max_iter=200, max_time=5, 
    multi_k=mnfp_instance.k_amount < 120
)
```

## Generating Figures

To save a figure of the instance, use the `print_mnfp_graph` function. Passing a solution will display it in the figures. 

```julia
print_mnfp_graph(mnfp_instance, solution=out.solution, file_path="all_sum.png", multi=false)
```

Be careful with the `multi` parameter. If set to `true` it will create *one file for each commodity* in the folder passed to `multi_folder` (`"mnfp_multi"` by default)

```julia
print_mnfp_graph(mnfp_instance, solution=out.solution, file_path="all_sum.png", multi=true, multi_folder="mnfp_multi")
```




### 