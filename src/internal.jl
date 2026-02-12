##region internal functionalities
function _plot_voxel(gridID::Array{Vector}, addRef::Bool=true)
    if isnothing(canvas)
        global canvas = plot([mesh3d(x=0, y=0, z=0)], blank_layout())
        display(canvas)
    else
        react!(canvas, [mesh3d(x=0, y=0, z=0)], blank_layout())
    end
    if addRef
        add_ref_axes!(canvas, [0, 0, 0], refLen[])
    end

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    nx = size(gridID, 1)
    ny = size(gridID, 2)
    nz = size(gridID, 3)

    xmin = voxel.start[1]
    ymin = voxel.start[2]
    zmin = voxel.start[3]

    id_list = sort(unique(gridID))
    filter!(x -> x != [], id_list)

    @all pts1 pts2 pts3 pts4 pts5 pts6 pts7 pts8 = fill(0.0, 3)
    @all r g b = 0.0
    for ind in eachindex(id_list)
        if idDict[id_list[ind][end]] != 0
            ptsArray = []
            for i in 1:nx, j in 1:ny, k in 1:nz
                if gridID[i, j, k] == id_list[ind]
                    pts1 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts2 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts3 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts4 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 1.5) * dz + zmin]
                    pts5 = [(i - 1.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts6 = [(i - 0.5) * dx + xmin, (j - 1.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts7 = [(i - 0.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]
                    pts8 = [(i - 1.5) * dx + xmin, (j - 0.5) * dy + ymin, (k - 0.5) * dz + zmin]
    
                    if i == 1 || gridID[i-1, j, k] == [] || idDict[gridID[i-1, j, k][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts4)
                        push!(ptsArray, pts8)
                        push!(ptsArray, pts5)
                    end
                    if i == nx || gridID[i+1, j, k] == [] || idDict[gridID[i+1, j, k][end]] == 0
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts6)
                    end
                    if j == 1 || gridID[i, j-1, k] == [] || idDict[gridID[i, j-1, k][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts6)
                        push!(ptsArray, pts5)
                    end
                    if j == ny || gridID[i, j+1, k] == [] || idDict[gridID[i, j+1, k][end]] == 0
                        push!(ptsArray, pts4)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts8)
                    end
                    if k == 1 || gridID[i, j, k-1] == [] || idDict[gridID[i, j, k-1][end]] == 0
                        push!(ptsArray, pts1)
                        push!(ptsArray, pts2)
                        push!(ptsArray, pts3)
                        push!(ptsArray, pts4)
                    end
                    if k == nz || gridID[i, j, k+1] == [] || idDict[gridID[i, j, k+1][end]] == 0
                        push!(ptsArray, pts5)
                        push!(ptsArray, pts6)
                        push!(ptsArray, pts7)
                        push!(ptsArray, pts8)
                    end
                end
            end
            if !haskey(colorDict, idDict[id_list[ind][end]])
                @all r g b = round(Int, rand() * 255)
                colorDict[idDict[id_list[ind][end]]] = "rgb($r, $g, $b)"
            end
            voxel_obj = polygons(ptsArray, 4, colorDict[idDict[id_list[ind][end]]])
    
            addtraces!(canvas, voxel_obj)
            sleep(0.1)
        end
    end
end

function _add_geom(geo::Geometry, gridID::Array{Vector})
    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    np = length(geo.pos)
    xrange = getindex.(geo.pos, 1)
    yrange = getindex.(geo.pos, 2)
    zrange = getindex.(geo.pos, 3)

    if isempty(gridID)
        @all ngx ngy ngz = 1
    else
        ngx = size(gridID, 1)
        ngy = size(gridID, 2)
        ngz = size(gridID, 3)
    end

    xmin = _round((minimum([minimum(xrange), voxel.start[1]]) - shift[]*dx/2)/dx)*dx + shift[]*dx/2
    ymin = _round((minimum([minimum(yrange), voxel.start[2]]) - shift[]*dy/2)/dy)*dy + shift[]*dy/2
    zmin = _round((minimum([minimum(zrange), voxel.start[3]]) - shift[]*dz/2)/dz)*dz + shift[]*dz/2

    xmax = _round((maximum([maximum(xrange), voxel.start[1] + (ngx - 1) * dx]) - shift[]*dx/2)/dx)*dx + shift[]*dx/2
    ymax = _round((maximum([maximum(yrange), voxel.start[2] + (ngy - 1) * dy]) - shift[]*dy/2)/dy)*dy + shift[]*dy/2
    zmax = _round((maximum([maximum(zrange), voxel.start[3] + (ngz - 1) * dz]) - shift[]*dz/2)/dz)*dz + shift[]*dz/2

    x = collect(xmin:dx:xmax)
    y = collect(ymin:dy:ymax)
    z = collect(zmin:dz:zmax)

    nx = round(Int, (xmax - xmin) / dx + 1)
    ny = round(Int, (ymax - ymin) / dy + 1)
    nz = round(Int, (zmax - zmin) / dz + 1)
    
    gridID_new = Array{Vector}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        gridID_new[i, j, k] = []
    end
    
    if !isempty(gridID)
        for i in 1:ngx, j in 1:ngy, k in 1:ngz
            gridID_new[findfirst(x .== voxel.start[1])+(i-1), findfirst(y .== voxel.start[2])+(j-1), findfirst(z .== voxel.start[3])+(k-1)] = gridID[i, j, k]
        end
    end

    for n in 1:np
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            if !(geo.ID in gridID_new[indx[i], indy[j], indz[k]])
                push!(gridID_new[indx[i], indy[j], indz[k]], geo.ID)
            end
        end
    end

    global gridID = gridID_new

    voxel.start[1] = xmin
    voxel.start[2] = ymin
    voxel.start[3] = zmin
    
    voxel.grid = export_grid()
    
    return nothing
end

function _add_geom(geoList::Vector{Geometry}, gridID::Array{Vector})
    for i in eachindex(geoList)
        _add_geom(geoList[i], gridID)
    end
end

function _del_geom(geo::Geometry, gridID::Array{Vector}, trim::Bool=true)
    np = length(geo.pos)

    dx = voxel.dl[1]
    dy = voxel.dl[2]
    dz = voxel.dl[3]

    ngx = size(gridID, 1)
    ngy = size(gridID, 2)
    ngz = size(gridID, 3)

    x = collect(voxel.start[1]:dx:(ngx-1)*dx+voxel.start[1])
    y = collect(voxel.start[2]:dy:(ngy-1)*dy+voxel.start[2])
    z = collect(voxel.start[3]:dz:(ngz-1)*dz+voxel.start[3])

    for n in 1:np
        indx = _find_nearest(x, geo.pos[n][1])
        indy = _find_nearest(y, geo.pos[n][2])
        indz = _find_nearest(z, geo.pos[n][3])
        for i in eachindex(indx), j in eachindex(indy), k in eachindex(indz)
            filter!(e -> e != geo.ID, gridID[indx[i], indy[j], indz[k]])
        end
    end

    if trim
        x1 = 1
        x2 = ngx
        y1 = 1
        y2 = ngy
        z1 = 1
        z2 = ngz

        for i in 1:ngx
            if unique(gridID[i, :, :]) != []
                x1 = i
                break
            end
        end
        for i in ngx:1
            if unique(gridID[i, :, :]) != []
                x2 = i
                break
            end
        end
        for i in 1:ngy
            if unique(gridID[:, i, :]) != []
                y1 = i
                break
            end
        end
        for i in ngy:1
            if unique(gridID[:, i, :]) != []
                y2 = i
                break
            end
        end
        for i in 1:ngz
            if unique(gridID[:, :, i]) != []
                z1 = i
                break
            end
        end
        for i in ngz:1
            if unique(gridID[:, :, i]) != []
                z2 = i
                break
            end
        end

        gridID = gridID[x1:x2, y1:y2, z1:z2]
        voxel.start[1] = voxel.start[1] + (x1 - 1) * dx
        voxel.start[2] = voxel.start[2] + (y1 - 1) * dy
        voxel.start[3] = voxel.start[3] + (z1 - 1) * dz
        
        voxel.grid = export_grid()
        
        return nothing
    end
end

function _del_geom(geoList::Vector{Geometry}, gridID::Array{Vector}, trim::Bool=true)
    for i in eachindex(geoList)
        _del_geom(geoList[i], gridID, trim)
    end
end

function _reset_gridID()
    empty!(idDict)
    idCount[] = 0
    grid_ind = sort(unique(voxel.grid))
    filter!(x -> x != 0, grid_ind)
    
    for ind in grid_ind
        idCount[] += 1
        idDict[idCount[]] = ind
    end
    nx = size(voxel.grid, 1)
    ny = size(voxel.grid, 2)
    nz = size(voxel.grid, 3)
    global gridID = Array{Vector}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        gridID[i, j, k] = []
        if voxel.grid[i, j, k] != 0
            ind = findfirst(x -> x .== voxel.grid[i, j, k], collect(values(idDict)))
            push!(gridID[i, j, k], collect(keys(idDict))[ind])
        end
    end
end
function _find_nearest(ary::Array{<:Number}, ele::Number)
    tmp = abs.(ary .- ele)
    val, ind = findmin(tmp)
    ind = findall(tmp .== val)
    return ind
end

function _round(num::Number)
    if num - floor(num) == 0.5
        if num < 0
            return floor(num)     
        else
            return ceil(num)
        end
    else
        return round(num)
    end
end

struct _TriangleData
    v0::NTuple{3, Float64}
    e1::NTuple{3, Float64}
    e2::NTuple{3, Float64}
    minv::NTuple{3, Float64}
    maxv::NTuple{3, Float64}
end

_sub3(a::NTuple{3, Float64}, b::NTuple{3, Float64}) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
_dot3(a::NTuple{3, Float64}, b::NTuple{3, Float64}) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross3(a::NTuple{3, Float64}, b::NTuple{3, Float64}) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function _auto_grid_coords(vmin::Float64, vmax::Float64, n::Int)
    @assert n > 0
    if n == 1
        return [(vmin + vmax) / 2], 1.0
    end
    step = (vmax - vmin) / (n + 0.5)
    vals = collect(range(vmin + step / 2, vmax - step / 2, length=n))
    return vals, step
end

function _to_triangle_data(tris::Vector{NTuple{3, NTuple{3, Float64}}})
    data = Vector{_TriangleData}(undef, length(tris))
    for i in eachindex(tris)
        v0 = tris[i][1]
        v1 = tris[i][2]
        v2 = tris[i][3]
        e1 = _sub3(v1, v0)
        e2 = _sub3(v2, v0)
        minv = (min(v0[1], v1[1], v2[1]), min(v0[2], v1[2], v2[2]), min(v0[3], v1[3], v2[3]))
        maxv = (max(v0[1], v1[1], v2[1]), max(v0[2], v1[2], v2[2]), max(v0[3], v1[3], v2[3]))
        data[i] = _TriangleData(v0, e1, e2, minv, maxv)
    end
    return data
end

function _unique_count_sorted!(vals::Vector{Float64}, atol::Float64)
    isempty(vals) && return 0
    sort!(vals)
    count = 1
    prev = vals[1]
    for i in 2:length(vals)
        if abs(vals[i] - prev) > atol
            count += 1
            prev = vals[i]
        end
    end
    return count
end

function _ray_triangle_t(origin::NTuple{3, Float64}, dir::NTuple{3, Float64}, tri::_TriangleData, eps::Float64)
    pvec = _cross3(dir, tri.e2)
    det = _dot3(tri.e1, pvec)
    if abs(det) < eps
        return NaN
    end
    invdet = 1.0 / det
    tvec = _sub3(origin, tri.v0)
    u = _dot3(tvec, pvec) * invdet
    if u < -eps || u > 1 + eps
        return NaN
    end
    qvec = _cross3(tvec, tri.e1)
    v = _dot3(dir, qvec) * invdet
    if v < -eps || u + v > 1 + eps
        return NaN
    end
    t = _dot3(tri.e2, qvec) * invdet
    return t > eps ? t : NaN
end

function _point_inside_direction(point::NTuple{3, Float64}, tris::Vector{_TriangleData}, dirsym::Char; eps::Float64=1e-9)
    dir = dirsym == 'x' ? (1.0, 0.0, 0.0) : dirsym == 'y' ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0)
    ts = Float64[]
    sizehint!(ts, min(64, length(tris)))

    for tri in tris
        if dirsym == 'x'
            if point[2] < tri.minv[2] - eps || point[2] > tri.maxv[2] + eps || point[3] < tri.minv[3] - eps || point[3] > tri.maxv[3] + eps
                continue
            end
        elseif dirsym == 'y'
            if point[1] < tri.minv[1] - eps || point[1] > tri.maxv[1] + eps || point[3] < tri.minv[3] - eps || point[3] > tri.maxv[3] + eps
                continue
            end
        else
            if point[1] < tri.minv[1] - eps || point[1] > tri.maxv[1] + eps || point[2] < tri.minv[2] - eps || point[2] > tri.maxv[2] + eps
                continue
            end
        end

        t = _ray_triangle_t(point, dir, tri, eps)
        if isfinite(t)
            push!(ts, t)
        end
    end
    isodd(_unique_count_sorted!(ts, 1e-7))
end

function _parse_raydirections(raydirection::AbstractString)
    dirs = unique(collect(lowercase(raydirection)))
    filter!(d -> d in ('x', 'y', 'z'), dirs)
    @assert !isempty(dirs) "raydirection must include x, y, and/or z"
    return dirs
end

function _stl_format(fileName::AbstractString)
    sz = filesize(fileName)
    open(fileName, "r") do io
        head = read(io, min(80, sz))
        headtxt = lowercase(strip(String(Char.(head))))
        starts_solid = startswith(headtxt, "solid")
        binary_sized = sz >= 84 && (sz - 84) % 50 == 0
        if binary_sized && !starts_solid
            return :binary
        end
        if starts_solid
            seek(io, max(0, sz - 256))
            tail = read(io, min(256, sz))
            tailtxt = lowercase(String(Char.(tail)))
            return occursin("endsolid", tailtxt) ? :ascii : :binary
        end
        return binary_sized ? :binary : :ascii
    end
end

function _read_stl_ascii(fileName::AbstractString)
    triangles = NTuple{3, NTuple{3, Float64}}[]
    verts = NTuple{3, Float64}[]
    open(fileName, "r") do io
        for line in eachline(io)
            s = strip(line)
            if startswith(lowercase(s), "vertex")
                parts = split(s)
                if length(parts) >= 4
                    x = parse(Float64, parts[2])
                    y = parse(Float64, parts[3])
                    z = parse(Float64, parts[4])
                    push!(verts, (x, y, z))
                    if length(verts) == 3
                        push!(triangles, (verts[1], verts[2], verts[3]))
                        empty!(verts)
                    end
                end
            end
        end
    end
    return triangles
end

function _read_stl_binary(fileName::AbstractString)
    triangles = NTuple{3, NTuple{3, Float64}}[]
    open(fileName, "r") do io
        read(io, 80)
        nfacet = Int(read(io, UInt32))
        sizehint!(triangles, nfacet)
        for _ in 1:nfacet
            read(io, Float32)  # normal x
            read(io, Float32)  # normal y
            read(io, Float32)  # normal z
            v1 = (Float64(read(io, Float32)), Float64(read(io, Float32)), Float64(read(io, Float32)))
            v2 = (Float64(read(io, Float32)), Float64(read(io, Float32)), Float64(read(io, Float32)))
            v3 = (Float64(read(io, Float32)), Float64(read(io, Float32)), Float64(read(io, Float32)))
            read(io, UInt16)  # attribute byte count
            push!(triangles, (v1, v2, v3))
        end
    end
    return triangles
end

function _read_stl(fileName::AbstractString)
    fmt = _stl_format(fileName)
    tris = fmt == :binary ? _read_stl_binary(fileName) : _read_stl_ascii(fileName)
    @assert !isempty(tris) "No triangles found in STL file: $fileName"
    return tris
end
#endregion
