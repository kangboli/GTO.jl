using SpecialFunctions
using LinearAlgebra

abstract type AbstractGaussian end

struct CartesianGaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
end

"""
A gaussian is

∂ⁱ/∂xⁱ ∂ʲ/∂yʲ ∂ᵏ/∂zᵏ exp(-α(r - c)²)
"""
struct Gaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
end

i(g::Gaussian)::Int = g.i
j(g::Gaussian)::Int = g.j
k(g::Gaussian)::Int = g.k
α(g::Gaussian)::Float64 = g.α
c(g::Gaussian)::Vector{Float64} = g.c


function nuclear_potential(g::Gaussian, r::Vector{Float64})
    (2π / α(g)) R(i(g), j(g), k(g), 0, α(g), (r-c(g))...)
end

function R(i::Int, j::Int, k::Int, n::Int, α::Float64, a::Float64, b::Float64, c::Float64)
   T = α * (a^2 + b^2 + c^2)
   i == j == k == 0 && return (-2α)^n * F(n, T)
   rest = (α, a, b, c)
   i == j == 0 && return c * R(0, 0, k, n+1, rest...)  + k * R(0,0,k-1,n+1, rest...)
   i == 0 && return b * R(0, j, k, n+1, rest...) + j * R(0, j-1,k,n+1, rest...)
   return a * R(i,j,k,n+1, rest...) + i * R(i-1,j,k,n+1, rest...)
end


function F(n::Int, T::Float64)
    return 1/2 * T^(-1/2-n) * (gamma(1/2 + n) - gamma(1/2+n, T))
end


function d(n_a, n_b, N, α, pa, pb)
    rest = α, pa, pb
    N > n_a + n_b && return 0
    n_a != 0 && 
    return (1/(2 * α)) * d(n_a-1, n_b, N-1, rest...) + 
        pa * d(n_a-1,n_b,N,rest...) + (N+1) d(n_a-1,n_b,N+1,rest...)

    n_b != 0 && 
    return (1/(2 * α)) * d(n_a, n_b-1, N-1, rest...) + 
        pa * d(n_a,n_b-1,N,rest...) + (N+1) d(n_a,n_b-1,N+1,rest...)

    return 1
end

function gaussian_product(a::CartesianGaussian, b::CartesianGaussian)
    cp = c(a) * α(a) + c(b) * α(b)
    αp = α(a) + α(b)
    e = exp(-α(a) * α(b) / (α(a) + α(b)) * norm(c(a) - c(b))^2)
    Ni, Nj, Nk = i(a) + i(b), j(a) + j(b), k(a) + k(b)
    dis = [d(i(a), i(b), N, αp, cp-c(a), cp-c(b)) for N = 1:Ni]
    djs = [d(j(a), j(b), N, αp, cp-c(a), cp-c(b)) for N = 1:Nj]
    dks = [d(k(a), k(b), N, αp, cp-c(a), cp-c(b)) for N = 1:Nk]

    [((e *  dis[ni] * djs[nj] * dks[nk]),
        Gaussian(ni, nj, nk, αp, cp))
        for ni in 1:Ni, nj in 1:Nj, nk in 1:Nk]
end

function overlap_integral(g::Gaussian)
    (i(g) == 0 && j(g) == 0 && k(g) == 0) || return 0
    return (π / α(g))^(3/2)
end

function two_electron_integral(p::Gaussian, q::Gaussian)
    λ = 2 π^(5/2) * α(p)^(-1)  * α(q)^(-1) * (α(p) + α(q))^(-1/2)
    αc = α(p) * α(q) / (α(p) + α(q)) * norm(c(p) - c(q))^2
    pq = c(p) - c(q)
    λ * (-1)^(i(q) + j(q) + k(q)) * R(i(p) + i(q), j(p) + j(q), k(p) + k(q), αc..., pq...)
end


