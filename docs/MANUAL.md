# VoxelModel.jl Manual

## Overview

`VoxelModel.jl` is for building, editing, and visualizing voxelized 3D models.

Core ideas:
- The package keeps a global voxel workspace (`voxel`) with:
  - `grid::Array{Int,3}`: material/color index per voxel (`0` means empty)
  - `dl::Vector{Float64}`: voxel spacing `[dx, dy, dz]`
  - `start::Vector{Float64}`: coordinate of the first voxel center in each axis
- Geometry creation APIs (`create_*`) return `Geometry` objects that can be transformed (`trans!`, `rot!`) or removed (`clear_geom`).
- Rendering is handled with Plotly through `plot_voxel()`.

## Quick Start

```julia
using VoxelModel

reset_voxel()
reset_dl(1.0)

c = create_cube([0, 0, 0], 4, 1, "corner")
s = create_sphere([5, 5, 5], 3, 2)

plot_voxel()
```

## Coordinate and Grid Conventions

- `shift[]` controls center alignment of voxel centers.
  - `shift[] == true` (default): centers are offset by half a cell.
  - `shift[] == false`: centers lie on integer grid lines (for `dl = 1`).
- `reset_start(start)` snaps the start coordinate up to the nearest valid grid center.
- `reset_dl(...)` and `reset_shift(...)` reset the voxel workspace.

## Public Globals

### `colorDict`

```julia
colorDict::Dict{Int,String}
```

Maps material index to plot color. Example:

```julia
colorDict[4] = "pink"
```

If an index is missing, a random color is assigned when plotting.

### `canvas`

```julia
canvas
```

Current Plotly canvas object. Useful if you want to add custom traces.

## API Reference

### Workspace Management

#### `reset_voxel(; render=false)`

Reset the whole voxel workspace.

- Returns: `nothing`
- Keyword:
  - `render`: if `true`, immediately redraw

#### `export_voxel()`

Return a copy of current voxel state.

- Returns: `Voxels`

#### `save_voxel(fileName::String)`

Save voxel state to JLD.

- Returns: `nothing`

#### `load_voxel(fileName::String)`

Load voxel state from JLD (replaces current workspace).

- Returns: `nothing`

#### `plot_voxel(addRef::Bool=true)`

Render current voxel model.

- Argument:
  - `addRef`: draw reference axes if `true`
- Returns: `nothing`

#### `assign_voxel(grid::Array{Int,3}, dl::Vector{<:Real}=[1.0,1.0,1.0], start::Vector{<:Real}=[0,0,0])`
#### `assign_voxel(grid::Array{Int,3}, dl::Real, start::Vector{<:Real}=[shift[]*0.5, shift[]*0.5, shift[]*0.5])`

Assign external voxel data into workspace.

- `grid`: voxel index array
- `dl`: spacing per axis (or scalar for isotropic spacing)
- `start`: first voxel center coordinate
- Returns: `nothing`

#### `export_grid()`

Export the current effective index grid.

If geometries overlap, the latest-added geometry wins at each voxel.

- Returns: `Array{Int,3}`

### Grid Controls

#### `reset_ref(b::Bool, len::Real=refLen[])`

Enable/disable reference axes in plots and set axis length.

- Returns: `nothing`

#### `reset_shift(b::Bool)`

Set grid shift mode (`shift[]`) and reset workspace.

- Returns: `nothing`

#### `reset_dl(dl::Vector{<:Real})`
#### `reset_dl(dl::Real)`

Set voxel spacing and reset workspace.

- Returns: `nothing`

#### `reset_start(start::Vector{<:Real})`

Set workspace start coordinate (snapped to valid voxel centers).

- Returns: `nothing`

### Geometry Creation

All creation APIs support `; render=false`.

#### `create_cuboid(origin, dim, ind=1, mode="corner", fac=2; render=false)`

Create cuboid geometry.

- `origin`: 3-vector
- `dim`: `[lx, ly, lz]`
- `ind`: material/color index
- `mode`: `"corner"` or `"center"`
- `fac`: internal sampling factor relative to grid spacing
- Returns: `Geometry`

#### `create_cube(origin, dim, ind=1, mode="corner", fac=2; render=false)`

Create cube geometry.

- Returns: `Geometry`

#### `create_sphere(origin, radius, ind=1, fac=2; render=false)`

Create sphere geometry.

- Returns: `Geometry`

#### `create_ellipsoid(origin, par, ind=1, fac=2; render=false)`

Create ellipsoid geometry.

- `par`: semi-axes `[a, b, c]`
- Returns: `Geometry`

#### `create_cylinder(origin, radius, height, ind=1, fac=2; render=false)`

Create cylinder geometry (axis along +z from `origin`).

- Returns: `Geometry`

### Geometry Editing

#### `trans!(geo::Geometry, dl::Vector{<:Real}; render=false)`

Translate an existing geometry.

- Returns: `nothing`

#### `rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0]; render=false)`

Rotate geometry by `ang` degrees around `axis`.

- If `origin == [0]`, geometry centroid is used.
- Returns: `nothing`

#### `clear_geom(geo::Geometry; render=false)`
#### `clear_geom(geoList::Vector{Geometry}; render=false)`

Remove one or multiple geometries.

- Returns: `nothing`

### STL Voxelization

#### `voxelize_stl(fileName::String, gridN::Union{Int,NTuple{3,Int}}=100, ind::Int=1, raydirection::String="xyz"; render=false)`

Voxelize a watertight STL mesh (ASCII or binary) into current workspace.

- `fileName`: path to STL
- `gridN`: number of voxels per axis (`Int` or `(nx, ny, nz)`)
- `ind`: index assigned to inside voxels
- `raydirection`: any combination of `"x"`, `"y"`, `"z"`
- The workspace is reset and replaced by the voxelized result.
- Returns: `Voxels`

Example:

```julia
using VoxelModel

stl_path = joinpath(@__DIR__, "..", "examples", "sample.stl")
voxelize_stl(stl_path, (100, 100, 100), 1, "xyz")
plot_voxel(true)
```

See also `examples/ex_stl_voxelize.jl`.

## Notes and Tips

- For faster interactive work, build with `render=false`, then call `plot_voxel()` once.
- Use lower `fac` and/or coarser grid for speed; increase for better shape fidelity.
- Index `0` is empty space. Using `ind=0` effectively carves/removes regions.
- STL voxelization expects a closed (watertight) mesh for reliable inside/outside classification.
