neighbors2d(I...) = neighbors2d(I)
function neighbors2d(I::Tuple{Int,Int})
    neigh = [Vector{Int}(undef, 4) for i in 1:prod(I)]
    Ly = I[1]
    Lx = I[2]
    for ix in 1:Lx
        for iy in 1:Ly        
            ctr = iy + Ly*(ix-1)
            N = mod1(iy - 1, Ly), ix 
            E = iy, mod1(ix + 1, Lx)
            S = mod1(iy + 1, Ly), ix
            W = iy, mod1(ix - 1, Lx)
            neigh[ctr][1] = N[1] + Ly*(N[2]-1)
            neigh[ctr][2] = W[1] + Ly*(W[2]-1)
            neigh[ctr][3] = S[1] + Ly*(S[2]-1)
            neigh[ctr][4] = E[1] + Ly*(E[2]-1)
        end
    end
    neigh
end
