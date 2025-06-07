function mcmcslow(x::Lattice,niter::Int, β::Float64, h::Vector{Float64})
    n = prod(x.dims)
    spin = rand((-1,1),n)
    res = Matrix{Int}(undef,n,niter)
    for t in 1:niter
        one_mc_sweep!(x,spin,β,h)
        #println(t," ",sum(spin)/n)
        for i in x.site
            res[i,t] = spin[i]            
        end
    end
    res
end
function one_mc_sweep!(x,spin,β,h)
    idxperm = randperm(length(spin))
    for i in idxperm
        one_mc_step!(x,spin,i,β,h)
    end
end

function one_mc_step!(x,spin,site,β,h)
    Eold = energyslow(x,h,spin)
    spin[site] = -spin[site]
    Enew = energyslow(x,h,spin)
    ΔE = -Eold + Enew
    if ΔE <= 0 || rand() < exp( -β * ΔE )
        return
    else
        spin[site] = -spin[site]
    end
end

function energyslow(x,h,spin)
    e = 0.0
    @inbounds for i in eachindex(spin)
        sumj = 0 
        @simd for j in x.neig[i]
            sumj += spin[j]
        end
        e += (0.5*sumj + h[i])*spin[i]
    end
    return -e
end