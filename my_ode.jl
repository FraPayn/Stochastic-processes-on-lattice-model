



abstract type Odesolver end
struct Euler <: Odesolver end
struct Euler_cromer <: Odesolver end
struct Midpoint <: Odesolver end
struct Leapfrog <: Odesolver end

struct Odeproblem{T<:AbstractFloat}
    ϕ::Function
    x0::T
    v0::T
    Δt::T
    tspan::Tuple{T,T}
end

struct odesolution{T<:AbstractFloat}
    t::Vector{T}
    x::Vector{T}
    v::Vector{T}
end

function update(solver::Euler,ϕ,v,x,Δt)
    dv = v + ϕ(x)*Δt
    dx = x + v *Δt
    return dv,dx
end

function update(solver::Euler_cromer,ϕ,v,x,Δt)
    dv = v + ϕ(x)*Δt
    dx = x + dv * Δt
    return dv,dx
end

function update(solver::Midpoint,ϕ,v,x,Δt)
    dv = v + ϕ(x)*Δt
    dx = x + (v + dv)*Δt/2
    return dv,dx
end

function update(solver::Leapfrog,ϕ,v,x,Δt)
    dx = x + v * Δt
    dv = v + ϕ(dx)*Δt
    return dv,dx
end

function _solve(p::Odeproblem{T},s::T,v) where {T<:AbstractFloat}
    ϕ,x0,v0,Δt,tmin,tmax=p.ϕ,p.x0,p.v0,p.Δt,p.tspan...
    x=[x0]
    t=[tmin]
    v=[v0]
    print(Δt)
    if s == Leapfrog()
        v0 += ϕ(x0)*Δt/2 
    end
    for k in tmin+Δt:Δt:tmax
        dv, dx = update(s, ϕ, v[end], x[end], Δt)
        push!(x, dx)
        push!(v, dv)
        push!(t, k)
    end
    return t,x,v
end

function solve(solver::T,ϕ,x0,v0,tspan;Δt::Real=0.1, verbose::Bool=true) where T <: Odesolver

    problem = Odeproblem(ϕ,v0,x0,Δt,tspan)
    t,x,v = _solve(problem,solver,verbose)
    return odesolution(t,x,v)    
end

function energy(x::Vector,v::Vector,ω)
    E = 0.5*v.^2 + 0.5*ω*x.^2
    δ = ([E[i].-E[end] for i in 1:length(E)])./E[end]
    return δ
end