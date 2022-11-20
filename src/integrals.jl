using SpecialFunctions
using LinearAlgebra
using IterTools: product


export AbstractGaussian,
    CartesianGaussian,
    Gaussian,
    i,
    j,
    k,
    α,
    c,
    nuclear_potential,
    R_tensor,
    F,
    product_coefficients,
    gaussian_product,
    two_electron_integral,
    overlap_integral,
    hermite,
    Lambda,
    evaluate,
    gaussians,
    coefficients,
    SphericalGaussian,
    l,
    m


abstract type AbstractGaussian end

"""
A Cartesian gaussian is

    (x-c₁)ⁱ (y-c₁)ʲ (z-c₁)ᵏ exp(-α(r-c)^2)
"""

struct CartesianGaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
    cache_number::Int
end

function CartesianGaussian(i,j,k,α,c)
    return CartesianGaussian(i,j,k,α,c,0) 
end


"""
A gaussian is

  ∂ⁱ/∂c₁ⁱ ∂ʲ/∂c₂ʲ ∂ᵏ/∂c₃ᵏ exp(-α(r - c)²)
= Λᵢ(x, α) Λⱼ(y, α) Λₖ(z, α) exp(-α(r - c)²)
"""
struct Gaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
end

struct SphericalGaussian <: AbstractGaussian
    l::Int
    m::Int
    α::Float64
    c::Vector{Float64}
end

l(g::SphericalGaussian) = g.l
m(g::SphericalGaussian) = g.m

struct ContractedGaussian
    coefficients::Vector{Float64}
    gaussians::Vector{<:AbstractGaussian}
end

coefficients(g::ContractedGaussian) = g.coefficients
gaussians(g::ContractedGaussian) = g.gaussians


i(g::AbstractGaussian)::Int = g.i
j(g::AbstractGaussian)::Int = g.j
k(g::AbstractGaussian)::Int = g.k
α(g::AbstractGaussian)::Float64 = g.α
c(g::AbstractGaussian)::Vector{Float64} = g.c

function nuclear_potential(cg::ContractedGaussian, r::Vector{<:Real})
    sum(zip(coefficients(cg), gaussians(cg))) do (c, g)
        c * nuclear_potential(g, r)
    end
end


function nuclear_potential(g::Gaussian, r::Vector{<: Real})
    (2π / α(g)) * R_tensor(i(g), j(g), k(g), 0, α(g), (r-c(g))...)
end

function R_tensor(i::Int, j::Int, k::Int, n::Int, α::Float64, a::Float64, b::Float64, c::Float64)
   T = α * (a^2 + b^2 + c^2)
   i == j == k == 0 && return (-2α)^n * F(n, T)
   rest = (α, a, b, c)
   R(i::Int, j::Int, k::Int, n::Int) = R_tensor(i,j,k,n, rest...)
   i == j == 0 && return c * R(0, 0, k, n+1)  + k * R(0,0,k-1,n+1)
   i == 0 && return b * R(0, j, k, n+1) + j * R(0, j-1,k,n+1)
   return a * R(i,j,k,n+1) + i * R(i-1,j,k,n+1)
end


function F(n::Int, T::Float64)
    return 1/2 * T^(-1/2-n) * (gamma(1/2+n) - gamma(1/2+n, T))
end


function product_coefficients(n_a, n_b, N, α, pa, pb)
    # println("$(n_a), $(n_b), $(N)")
    rest = α, pa, pb
    d(n_a, n_b, N) = product_coefficients(n_a, n_b, N, rest...)
    0 <= N <= n_a + n_b || return 0
    n_a != 0 && return let n_a = n_a - 1 
        1/(2α) * d(n_a, n_b, N-1) + pa * d(n_a,n_b,N) + (N+1) * d(n_a,n_b,N+1)
    end

    n_b != 0 && return let n_b = n_b - 1
        1/(2α) * d(n_a, n_b, N-1) + pb * d(n_a,n_b,N) + (N+1) * d(n_a,n_b,N+1)
    end

    return 1
end

function gaussian_product(a::CartesianGaussian, b::CartesianGaussian)
    αp = α(a) + α(b)
    cp = (c(a) * α(a) + c(b) * α(b)) / αp
    e = exp(-α(a) * α(b) / αp * norm(c(a) - c(b))^2)
    Ni, Nj, Nk = i(a) + i(b), j(a) + j(b), k(a) + k(b)
    dis = [product_coefficients(i(a), i(b), N, αp, (cp-c(a))[1], (cp-c(b))[1]) for N = 0:Ni]
    djs = [product_coefficients(j(a), j(b), N, αp, (cp-c(a))[2], (cp-c(b))[2]) for N = 0:Nj]
    dks = [product_coefficients(k(a), k(b), N, αp, (cp-c(a))[3], (cp-c(b))[3]) for N = 0:Nk]

    coefficients_and_gaussians = [
        ((e * dis[ni + 1] * djs[nj + 1] * dks[nk + 1]), Gaussian(ni, nj, nk, αp, cp)) for
        ni = 0:Ni, nj = 0:Nj, nk = 0:Nk
        ] |> vec

    ContractedGaussian(first.(coefficients_and_gaussians), last.(coefficients_and_gaussians))
end

function overlap_integral(g::Gaussian)
    (i(g) == 0 && j(g) == 0 && k(g) == 0) || return 0
    return (π / α(g))^(3/2)
end

function two_electron_integral(p::Gaussian, q::Gaussian)
    λ = 2 * π^(5/2) * α(p)^(-1) * α(q)^(-1) * (α(p) + α(q))^(-1/2)
    αc = α(p) * α(q) / (α(p) + α(q)) * norm(c(p) - c(q))^2
    pq = c(p) - c(q)
    λ * (-1)^(i(q) + j(q) + k(q)) * R_tensor(i(p) + i(q), j(p) + j(q), k(p) + k(q), αc..., pq...)
end


function evaluate(g::CartesianGaussian, r::Vector{<:Real})
    x,y,z = r
    c1,c2,c3 = c(g)
    (x - c1)^i(g) * (y - c2)^j(g) * (z - c3)^k(g) * exp(-α(g) * norm(r-c(g))^2)
end

function evaluate(g::Gaussian, r::Vector{<:Real})
    x,y,z = r
    c1,c2,c3 = c(g)
    Lambda(i(g), x-c1, α(g)) * Lambda(j(g), y-c2, α(g)) * Lambda(k(g), z-c3, α(g)) * exp(-α(g) * norm(r - c(g))^2)
end

function evaluate(g::ContractedGaussian, r::Vector{<:Real})
    sum(zip(coefficients(g), gaussians(g))) do (c, g)
        c * evaluate(g, r)
    end
end

function hermite(N::Int, x::Number)
    N == 0 && return 1
    N == 1 && return 2x
    return 2x * hermite(N-1, x) - 2 * (N-1) * hermite(N-2, x)
end

function Lambda(N::Int, x, α) 
    α^(N/2) * hermite(N, α^(1/2) * x)
end

function double_factorial(n::Int)
    n <= 1 && return 1
    n * double_factorial(n - 2)
end

function normalization(g::CartesianGaussian)
    n(v) = (2α(g)/π)^(1/4) * (4α(g))^(v/2) * (double_factorial(2v-1))^(-1/2)
    prod(n, [i(g), j(g), k(g)])
end

function normalization(g::SphericalGaussian)
    n = l(g)
    (factorial(2n + 2) * π^(1/2) /
     (2^(2n+3) * factorial(n+1) * α^(n+3/2)))^(-1/2)
end

#= function spherical_to_cartesian(g::SphericalGaussian)
    f = factorial
    map(filter(ls->sum(ls) == l(g), collect(product(0:l(g), 0:l(g), 0:l(g))))) do (l_x, l_y, l_z)
        j = floor((l_x + l_y - abs(m))/2)
        sqrt(f(2l_x) * f(2l_y) * f(2l_z) * f(l(g)) * f(l(g) - abs(m(g))) / 
        f(2l) * f(l_x) * f(l_y) * f(l_z) * f(l(g) + abs(m(g)))) * 
        1/(2^(l(g)) * f(l(g))) * sum(0:floor((l(g)-abs(m(g)))/2)) do i
            binomial(l(g), i) * binomial(i, j)
        end
    end
end
 =#
