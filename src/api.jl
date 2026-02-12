#region API
"""
    reset_voxel(;render=false)

    Reset the full voxel voxel. 
    
    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function reset_voxel(;render=false)
    global voxel = Voxels()
    global gridID = []
    idCount[] = 0
    empty!(idDict)
    if render
        _plot_voxel(gridID, refAxis[])
    end
    return nothing
end

"""
    export_voxel()

    Exports a copy of the current voxel.
    
    # Returns
    - `voxelCopy::Voxels`: The copy of the current voxel space.
"""
function export_voxel()
    voxelCopy =  Voxels(voxel.grid, voxel.dl, voxel.start)
    return voxelCopy
end

"""
    save_voxel(fileName::String)

    save voxel in JLD format. 
    
    # Arguments
    - `fileName::String`: file name of the JLD file
"""
function save_voxel(fileName::String)
    save(fileName, "voxel", voxel)
    return nothing
end

"""
    load_voxel(fileName::String)

    load voxel in JLD format. This will reset the current voxel space.
    
    # Arguments
    - `fileName::String`: file name of the JLD file
"""
function load_voxel(fileName::String)
    global voxel = load(fileName, "voxel")
    _reset_gridID()
    
    return nothing
end


"""
    plot_voxel(addRef::Bool=true)

    Plots the voxel voxel. If `addRef=false`, the reference axes will not be added. Call this function if the plot window is accidentally closed.
    
    # Arguments
    - `addRef::Bool=true`: Boolean value to specify whether to add the reference axes to the plot.
"""
function plot_voxel(addRef::Bool=true)
    
    global canvas = plot([mesh3d(x=0, y=0, z=0)], blank_layout())
    display(canvas)
    
    _plot_voxel(gridID, addRef)
end

"""
assign_voxel(grid::Array{Int, 3}, dl::Vector{<:Real}=[1.0, 1.0, 1.0], start::Vector{<:Real}=[0, 0, 0])

    Assign grid to voxel space.
    
    # Arguments
    - `grid::Array{Int, 3}`: Interger grid array.
    - `dl::Vector{<:Real}=[1.0, 1.0, 1.0]`: grid spacings.
    - `start::Vector{<:Real}=[shift[] * 0.5, shift[] * 0.5, shift[] * 0.5]`: start point.
"""
function assign_voxel(grid::Array{Int, 3}, dl::Vector{<:Real}=[1.0, 1.0, 1.0], start::Vector{<:Real}=[0, 0, 0])
    reset_dl(dl)
    reset_start(start)
    global voxel.grid = grid
    _reset_gridID()
    return nothing
end
function assign_voxel(grid::Array{Int, 3}, dl::Real, start::Vector{<:Real}=[shift[] * 0.5, shift[] * 0.5, shift[] * 0.5])
    assign_voxel(grid, [dl, dl, dl], start)
    return nothing
end

"""
    reset_ref(b::Bool, len::Real=refLen[])

    Toggles the display of the reference axes at the origin. The default state is `true` (axes visible).
    
    # Arguments
    - `b::Bool`: Boolean value to set the visibility of the reference axes.
    - `len::Float64=refLen[]`: reference length of the axes. The default is the minimum of the grid spacings.
"""
function reset_ref(b::Bool, len::Real=refLen[])
    refAxis[] = b
    refLen[] = Float64(len)
    return nothing
end

"""
    reset_shift(b::Bool)

    Sets the `shift[]` parameter to the specified boolean value `b`. 
    !! Since the gird space is changed, the voxel space will be reseted accordingly. 
    Use this function before adding geomtries to the voxel space.
    
    # Arguments
    - `b::Bool`: Boolean value to set the `shift[]` parameter.
"""
function reset_shift(b::Bool)
    shift[] = b
    reset_voxel()
    return nothing
end

"""
    reset_dl(dl::Vector{<:Real})

    Updates the grid spacing to the specified vector `dl`.
    !! Since the gird space is changed, the voxel space will be reseted accordingly. 
    Use this function before adding geomtries to the voxel space.
    
    # Arguments
    - `dl::Vector{<:Real}`: A vector containing the new grid spacing values.
"""
function reset_dl(dl::Vector{<:Real})
    @assert length(dl) == 3
    @assert all(>(0), dl)
    reset_voxel()
    voxel.dl = dl
    reset_start([0.0, 0.0, 0.0])
    reset_ref(true, minimum(dl))
    return nothing
end
function reset_dl(dl::Real)
    reset_dl([dl, dl, dl])
    return nothing
end

"""
    reset_start(start::Vector{<:Real})

    reset the start point of the voxel space. 
    !! Note that th spart point is ceiled to the neasrest grid center.
    
    # Arguments
    - `start::Vector{<:Real}`: A vector indicating the spart point.
"""
function reset_start(start::Vector{<:Real})
    @assert length(start) == 3
    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]
    
    xmin = ceil((start[1] - shift[]*dx/2)/dx)*dx + shift[]*dx/2
    ymin = ceil((start[2] - shift[]*dy/2)/dy)*dy + shift[]*dy/2
    zmin = ceil((start[3] - shift[]*dz/2)/dz)*dz + shift[]*dz/2
    voxel.start = [xmin, ymin, zmin]
    return nothing
end

"""
    create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode::String="corner", fac::Real=2; render=false)

    Creates a cuboid with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the cuboid.
    - `dim::Vector{<:Real}`: The dimensions of the cuboid.
    - `ind::Int=1`: The color index of the cuboid.
    - `mode::String="corner"`: The mode specifying the cuboid's origin ("corner" or "center").
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
    
    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function create_cuboid(origin::Vector{<:Real}, dim::Vector{<:Real}, ind::Int=1, mode="corner", fac::Real=2; render=false)
    @assert length(origin) == 3
    @assert length(dim) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    sx = Int(_round((dim[1]-2*dx/fac) / (dx/fac))) + 1
    sy = Int(_round((dim[2]-2*dy/fac) / (dy/fac))) + 1
    sz = Int(_round((dim[3]-2*dz/fac) / (dz/fac))) + 1

    pos = []
    if mode == "center"
        xs = origin[1] - dim[1]/2 + dx / fac
        ys = origin[2] - dim[2]/2 + dy / fac
        zs = origin[3] - dim[3]/2 + dz / fac
    else
        xs = origin[1] + dx / fac
        ys = origin[2] + dy / fac
        zs = origin[3] + dz / fac
    end

    for i in 1:sx, j in 1:sy, k in 1:sz
        push!(pos, [xs + (i - 1) * dx/fac, ys + (j - 1) * dy/fac, zs + (k - 1) * dz/fac])
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    if render
        _plot_voxel(gridID, refAxis[])
    end
    return geo
end

"""
    create_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode::String="corner", fac::Real=2; render=false)

    Creates a cube with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the cube.
    - `dim::Real`: The side length of the cube.
    - `ind::Int=1`: The color index of the cube.
    - `mode::String="corner"`: The mode specifying the cube's origin ("corner" or "center").
    - `fac::Real=2`: The interior densified factor according to the grid spacing.
    
    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function create_cube(origin::Vector{<:Real}, dim::Real, ind::Int=1, mode="corner", fac::Real=2; render=false)
    @assert dim > 0
    return create_cuboid(origin, [dim, dim, dim], ind, mode, fac; render=render)
end

"""
    create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2; render=false)

    Creates a sphere with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the sphere.
    - `radius::Real`: The radius of the sphere.
    - `ind::Int=1`: The color index of the sphere.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
    """
function create_sphere(origin::Vector{<:Real}, radius::Real, ind::Int=1, fac::Real=2; render=false)
    @assert length(origin) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    sr = maximum([ceil(Int, radius / (dx/fac)), ceil(Int, radius / (dy/fac)), ceil(Int, radius / (dz/fac))])

    pos = []
    for i in -sr:sr, j in -sr:sr, k in -sr:sr

        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1])./(radius-dx/fac))^2 + ((y1 - origin[2])./(radius-dy/fac))^2 + ((z1 - origin[3])./(radius-dz/fac))^2 < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)

    if render
        _plot_voxel(gridID, refAxis[])
    end
    return geo
end

"""
    create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2; render=false)

    Creates an ellipsoid with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The origin point of the ellipsoid.
    - `par::Vector{<:Real}`: The lengths of the semi-axes.
    - `ind::Int=1`: The color index of the ellipsoid.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
    """
function create_ellipsoid(origin::Vector{<:Real}, par::Vector{<:Real}, ind::Int=1, fac::Real=2; render=false)
    @assert length(origin) == 3
    @assert length(par) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    sa = ceil(Int, par[1] / (dx/fac))
    sb = ceil(Int, par[2] / (dy/fac))
    sc = ceil(Int, par[3] / (dz/fac))

    pos = []
    for i in -sa:sa, j in -sb:sb, k in -sc:sc
        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1]) / (par[1]-dx/fac))^2 + ((y1 - origin[2]) / (par[2]-dy/fac))^2 + ((z1 - origin[3]) / (par[3]-dz/fac))^2 < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)
    if render
        _plot_voxel(gridID, refAxis[])
    end
    return geo
end

"""
    create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2; render=false)

    Creates a cylinder with the specified parameters.
    
    # Arguments
    - `origin::Vector{<:Real}`: The base origin point of the cylinder.
    - `radius::Real`: The radius of the cylinder.
    - `height::Real`: The height of the cylinder.
    - `ind::Int=1`: The color index of the cylinder.
    - `fac::Real=2`: The interior densified factor according to the grid spacing.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
    """
function create_cylinder(origin::Vector{<:Real}, radius::Real, height::Real, ind::Int=1, fac::Real=2; render=false)
    @assert length(origin) == 3
    @assert fac > 0

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    sr = maximum([ceil(Int, radius / (dx/fac)), ceil(Int, radius / (dy/fac))])
    sz = Int(_round((height-2*dz/fac) / (dz/fac))) + 1
    pos = []
    for i in -sr:sr, j in -sr:sr, k in 1:sz
        x1 = origin[1] + i * (dx/fac) 
        y1 = origin[2] + j * (dy/fac) 
        z1 = origin[3] + k * (dz/fac) 

        if ((x1 - origin[1])./(radius-dx/fac))^2 + ((y1 - origin[2])./(radius-dy/fac))^2  < 1
            push!(pos, [x1, y1, z1])
        end
    end

    idCount[] += 1
    geo = Geometry(pos, ind, idCount[])
    idDict[idCount[]] = ind
    _add_geom(geo, gridID)
    if render
        _plot_voxel(gridID, refAxis[])
    end
    return geo
end

"""
    voxelize_stl(fileName::String, gridN::Union{Int, NTuple{3, Int}}=100, ind::Int=1, raydirection::String="xyz"; render=false)

    Voxelizes a watertight STL mesh and assigns the result to the current voxel space.

    # Arguments
    - `fileName::String`: path to STL file (ASCII or binary).
    - `gridN::Union{Int, NTuple{3, Int}}=100`: voxel counts in `(x, y, z)`. If integer, same count is used for all dimensions.
    - `ind::Int=1`: color/material index for occupied voxels.
    - `raydirection::String="xyz"`: ray-casting direction(s), any combination of `"x"`, `"y"`, `"z"`.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.

    # Returns
    - `Voxels`: copy of the voxelized result.
"""
function voxelize_stl(fileName::String, gridN::Union{Int, NTuple{3, Int}}=100, ind::Int=1, raydirection::String="xyz"; render=false)
    @assert ind >= 0
    nx, ny, nz = gridN isa Int ? (gridN, gridN, gridN) : gridN
    @assert nx > 0 && ny > 0 && nz > 0

    tris = _to_triangle_data(_read_stl(fileName))
    xmin = minimum(getfield.(tris, :minv) .|> x -> x[1])
    ymin = minimum(getfield.(tris, :minv) .|> x -> x[2])
    zmin = minimum(getfield.(tris, :minv) .|> x -> x[3])
    xmax = maximum(getfield.(tris, :maxv) .|> x -> x[1])
    ymax = maximum(getfield.(tris, :maxv) .|> x -> x[2])
    zmax = maximum(getfield.(tris, :maxv) .|> x -> x[3])

    xs, dx = _auto_grid_coords(xmin, xmax, nx)
    ys, dy = _auto_grid_coords(ymin, ymax, ny)
    zs, dz = _auto_grid_coords(zmin, zmax, nz)
    dirs = _parse_raydirections(raydirection)

    inside_votes = zeros(UInt8, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        point = (xs[ix], ys[iy], zs[iz])
        for d in dirs
            if _point_inside_direction(point, tris, d)
                inside_votes[ix, iy, iz] += 1
            end
        end
    end

    threshold = ceil(Int, length(dirs) / 2)
    grid = zeros(Int, nx, ny, nz)
    for i in eachindex(grid)
        if inside_votes[i] >= threshold
            grid[i] = ind
        end
    end

    reset_voxel()
    voxel.grid = grid
    voxel.dl = [dx, dy, dz]
    voxel.start = [xs[1], ys[1], zs[1]]
    _reset_gridID()

    if render
        _plot_voxel(gridID, refAxis[])
    end
    return export_voxel()
end

"""
    trans!(geo::Geometry, dl::Vector{<:Real}; render=false)

    Translates the geometry by the specified vector `dl`.
    
    # Arguments
    - `geo::Geometry`: The geometry to be translated.
    - `dl::Vector{<:Real}`: The translation vector.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
    """
function trans!(geo::Geometry, dl::Vector{<:Real}; render=false)
    @assert length(dl) == 3

    _del_geom(geo, gridID)

    for n in eachindex(geo.pos)
        geo.pos[n][1] += dl[1]
        geo.pos[n][2] += dl[2]
        geo.pos[n][1] += dl[3]
    end

    _add_geom(geo, gridID)

    if render
        _plot_voxel(gridID, refAxis[])
    end
    return nothing
end

"""
    rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0]; render=false)

    Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.
    
    # Arguments
    - `geo::Geometry`: The geometry to be rotated.
    - `ang::Real`: The rotation angle.
    - `axis::Vector{<:Real}`: The rotation axis.
    - `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function rot!(geo::Geometry, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0]; render=false)
    @assert length(axis) == 3

    _del_geom(geo, gridID)

    axis = axis ./ norm(axis)
    vrot = similar(geo.pos)
    
    if origin == [0] # rotation center set at the geometry center
        origin = sum(geo.pos) ./ length(geo.pos)
    else
        @assert length(origin) == 3
    end

    for n in eachindex(vrot)
        v = (geo.pos[n] .- origin)
        vrot[n] = cosd(ang) * v + sind(ang) * cross(axis, v) + (1-cosd(ang)) * dot(axis, v) * axis
        geo.pos[n] = vrot[n] .+ origin
    end

    _add_geom(geo, gridID)

    if render
        _plot_voxel(gridID, refAxis[])
    end
    return nothing
end

"""
    clear_geom(geo::Geometry; render=false)

    Removes the specified geometry from the voxel voxel.
    
    # Arguments
    - `geo::Geometry`: The geometry to be removed.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function clear_geom(geo::Geometry; render=false)
    _del_geom(geo, gridID)
    geo = nothing
    if render
        _plot_voxel(gridID, refAxis[])
    end
    return nothing
end

"""
    clear_geom(geoList::Vector{Geometry}; render=false)

    Removes the specified list of geometries from the voxel voxel.
    
    # Arguments
    - `geoList::Vector{Geometry}`: The list of geometries to be removed.

    # Keywords
    - `render=false`: real-time rendering for creation/operation.
"""
function clear_geom(geoList::Vector{Geometry}; render=false)
    for i in eachindex(geoList)
        clear_geom(geoList[i], render=render)
        geoList[i] = nothing
    end
    return nothing
end

"""
    export_grid()

    Exports the grid array filled with color indexes. Note that when geometries overlap, the index of the last-added geometry is used.
    
    # Returns
    - `Array{Int}`: The grid array with color indexes.
"""
function export_grid()
    grid = zeros(Int, size(gridID))
    for i in eachindex(grid)
        if gridID[i] == []
            grid[i] = 0
        else
            grid[i] = idDict[gridID[i][end]]
        end
    end
    return grid
end


#endregion
