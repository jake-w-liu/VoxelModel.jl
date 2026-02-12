using VoxelModel

# Voxelize the STL mesh in this folder and visualize it.
stl_path = joinpath(@__DIR__, "sample.stl")

# Optional: ensure a clean scene.
reset_voxel()

# Voxelize using a 100x100x100 grid, fill index 1, robust xyz ray casting.
voxelize_stl(stl_path, (100, 100, 100), 1, "xyz")

# Show the voxelized model.
plot_voxel(true)
