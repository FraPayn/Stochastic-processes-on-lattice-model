function cartesian(ix,iy,Ny)
    return iy + Ny*(ix-1)
end

lin2cart(i, dims::NTuple) = CartesianIndices(dims)[i].I

function neighbors2d(I::Tuple{Int,Int})
    Nx=I[1]
    Ny=I[2]
    V = Vector{Vector{Int64}}()  
    E=N=S=W=0
    for i in 1:Nx
        for j in 1:Ny
            j == 1 ? N = cartesian(i,Ny,Ny) : N = cartesian(i,j-1,Ny)
            i == 1 ? W = cartesian(Nx,j,Ny) : W = cartesian(i-1,j,Ny) 
            j == Ny ?  S = cartesian(i,1,Ny) : S = cartesian(i,j+1,Ny)
            i == Nx ?  E = cartesian(1,j,Ny) : E = cartesian(i+1,j,Ny)
            push!(V,[N,W,S,E])
        end
    end
    return V 
end

function randomwalk1(dims::Tuple{Int,Int};x0::Tuple{Int,Int} = (1,1),titer::Int=1000)
    Nx=dims[1]
    Ny=dims[2]
    x_0 = cartesian(x0[1],x0[2],Ny)
    indices = Vector{Int64}()
    V = neighbors2d(dims)
    next_index = x_0
    push!(indices,x_0)
    for i in 1:titer
        s = rand()
        if 0.25 <= s <= 0.50
            next_index = V[next_index][1]
        end
        if 0.50 <= s <= 0.75
            next_index = V[next_index][4] 
        end
        if 0 <= s <= 0.25
            next_index = V[next_index][2] 
        end
        if 0.75 <= s <= 1
            next_index = V[next_index][3] 
        end
        push!(indices, next_index)

    end
    return indices
end

struct Lattice{D} 
    dims::NTuple{D,Int} # the tuple containing the dimensions 
    site::UnitRange{Int} # A range that tells me the sites 1 ... prod(dims)
    neig::Vector{Vector{Int}} # neigh[i] returns the vector of the neighbors of site i
end

function randomwalk(L::Lattice,nsteps::Int)
    d = length(L.dims)
    traj = zeros(Int, nsteps+1)
    traj[1] = cart2lin(L.dims,L.dims .>> 1 .+ 1) #lattice centre initialization saved as first trajectory point / .>> = *0.5
    for i in 2:nsteps+1
        traj[i] = L.neig[traj[i-1]][rand(1:2d)]
    end
    traj
end



(Lattice(dims::NTuple{D,Int}) where D) = Lattice(dims,1:prod(dims),neighbors(dims)) #constructor
Lattice(dims...) = Lattice(dims) # now I can call Lattice(3,5) instead of the default Lattice((3,5))

function show(io::IO, x::Lattice{D}) where D  # overload of the default print of the struct on stdout
    print(io,"$D-D lattice dims = $(x.dims[1])")
    for i in 2:length(x.dims)
        print(io,"×",x.dims[i])
    end
    print(io," PBC")
end


neighbors(I...) = neighbors(I)



function neighbors(I::NTuple{N,Int}) where N
    neigh = Vector{Vector{Int}}(undef,prod(I))
    vscra = zeros(Int,N)
    pm = (-1,1)
    site = 0
    @inbounds for i in CartesianIndices(I)
        cnt = 0
        site += 1
        neigh[site] = Vector{Int}(undef,2N)
        for pm1 in pm
            for k1=1:N
                cnt += 1
                for k2=1:N                                
                    vscra[k2] = mod1(i.I[k2] + pm1 * (k1 == k2), I[k2])
                end
                neigh[site][cnt] = cart2lin(I,vscra)
            end
        end
    end   
    neigh 
end

# utility to compute Linear index from a Cartesian representation
function cart2lin(I::NTuple{N,Int},v) where N 
    num = v[1]
    pdim  = 1
    @inbounds for i=1:N-1
        pdim *= I[i]
        num += pdim * (v[i+1]-1)
    end
    num
end

function trajdiffusion(vtraj::Vector{Vector{Int}}, L::Lattice)
    ntraj = length(vtraj)
    ave = zeros(length(vtraj[1]))
    va = similar(ave)
    for t in 1:length(vtraj[1])
        media = 0.0
        mom2 = 0.0
        for i in 1:ntraj
            p0 = lin2cart(vtraj[i][1],L.dims)
            pt = lin2cart(vtraj[i][t],L.dims)
            dt = sqrt(sum((p0 .- pt).^2))
            media += dt
            mom2 += dt^2
        end
        ave[t] = media/ntraj
        va[t] = (mom2-media*media/ntraj)/ntraj
    end
    ave,va
end

function trajRMSD(vtraj::Vector{Vector{Int}}, L::Lattice)
    ntraj = length(vtraj)
    ave = zeros(length(vtraj[1]))
    for t in 1:length(vtraj[1])
        media = 0.0
        for i in 1:ntraj
            p0 = lin2cart(vtraj[i][1],L.dims)
            pt = lin2cart(vtraj[i][t],L.dims)
            dt = sum((p0 .- pt).^2)
            media += dt
        end
        ave[t] = sqrt(media/ntraj)
    end
    ave
end
