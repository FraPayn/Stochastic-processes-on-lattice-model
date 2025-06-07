import Base.show # import Base.show for redefining the stdout printing of the struct
struct Lattice{D} 
    dims::NTuple{D,Int} # the tuple containing the dimensions 
    site::UnitRange{Int} # A range that tells me the sites 1 ... prod(dims)
    neig::Vector{Vector{Int}} # neigh[i] returns the vector of the neighbors of site i
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

# utility to translate Linear into Cartesian indices
lin2cart(i, dims::NTuple) = CartesianIndices(dims)[i].I

nothing
