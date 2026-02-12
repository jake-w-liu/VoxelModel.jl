# VoxelModel.jl

[![Build Status](https://github.com/akjake616/VoxelModel.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/akjake616/VoxelModel.jl/actions/workflows/CI.yml)

`VoxelModel.jl` is a Julia module for creating, manipulating, and visualizing 3D voxel geometries. It provides functionalities to create geometric shapes, add or delete them from a voxel grid, and visualize the grid using [`PlotlyJS.jl`](https://github.com/JuliaPlots/PlotlyJS.jl) and [`PlotlyGeometries.jl`](https://github.com/akjake616/PlotlyGeometries.jl). This project is dedicated to (hopefully) easy model creation in FDTD/PSTD simulations.

<p align="center">
  <img alt="VoxelModel.jl" src="./media/illus.png" width="50%" height="auto" />
</p>

## Installation

To install `VoxelModel.jl`, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.add("VoxelModel")
```

## Learn by Examples

Please go to the examples folder for quick understanding of how to use the module:


- `ex_basics.jl`: This example code illustrates the basic functionalities of the module. 

- `ex_rotation.jl`: See how rotation works.

- `ex_advanced.jl`: Create more complicated geometries. This example demonstrates the creation of a sphere hollowd by cylinders. The voxel space is saved as a JLD file.

- `ex_load.jl`: Load voxel space from JLD file and perform further editing.

- `ex_assign.jl`: Assign grid array to the voxel space.

- `ex_stl_voxelize.jl`: Voxelize an STL mesh (`examples/sample.stl`) and visualize the result.


## Usage

Please refer to the [user manual](./docs/MANUAL.md) in the docs folder.

### STL Voxelization

`VoxelModel.jl` supports voxelization of watertight ASCII or binary STL files:

```julia
using VoxelModel

stl_path = joinpath(@__DIR__, "examples", "sample.stl")
voxelize_stl(stl_path, (100, 100, 100), 1, "xyz")
plot_voxel(true)
```



