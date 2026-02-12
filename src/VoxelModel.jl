module VoxelModel

export reset_voxel, plot_voxel, export_voxel, save_voxel, load_voxel, assign_voxel
export reset_shift, reset_dl, reset_start, reset_ref
export create_cuboid, create_cube, create_sphere, create_ellipsoid, create_cylinder
export voxelize_stl
export trans!, rot!, export_grid
export clear_geom
export colorDict, canvas

using PlotlySupply
using PlotlyGeometries
using LinearAlgebra
using Combinatorics
using BatchAssign
using JLD

#region structs
mutable struct Geometry
    pos::Vector{Vector{Float64}}
    index::Int
    const ID::Int
end

Base.@kwdef mutable struct Voxels
    grid::Array{Int, 3} = zeros(1, 1, 1)
    dl::Vector{Float64} = [1.0, 1.0, 1.0]
    start::Vector{Float64} = [shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]
end
#endregion

#region constrols
const shift = Ref(true)
const refAxis = Ref(true)
const refLen = Ref{Float64}(1)
const idCount = Ref{Int}(0)
const idDict = Dict{Int, Int}()
const colorDict = Dict{Int, String}()

global voxel = Voxels()
global gridID::Array{Vector} = []
global canvas = nothing
#endregion

include("api.jl")
include("internal.jl")

end
