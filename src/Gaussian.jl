using Accessors

export AbstractGaussian,
    i,
    j,
    k,
    α,
    c,
    gaussians,
    coefficients,
    CartesianGaussian,
    HermiteGaussian,
    ContractedGaussian,
    contract,
    set_i,
    set_j,
    set_k,
    set_c,
    set_α,
    CCG,
    CHG,
    Gaussian,
    Conj,
    term


abstract type AbstractGaussian end

"""
A Cartesian (primitive) Gaussian is

    (x-c₁)ⁱ (y-c₂)ʲ (z-c₃)ᵏ exp(-α(r-c)^2)
"""
struct CartesianGaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
    gaussian_cache_number::Int
end

CartesianGaussian(i, j, k, α, c) = CartesianGaussian(i, j, k, α, c, -1)
cache_number(g::CartesianGaussian) = g.gaussian_cache_number

function cache!(c::Cache, g::CartesianGaussian, atom::AbstractAtom)
    g = Accessors.@set g.gaussian_cache_number = genreate_cartesian_cache_number!(c)
    push!(c.cartesian_normalization_factors, normalization(g))
    push!(c.gaussian_atom, atom)
    return g
end

function Base.hash(g::CartesianGaussian, h::UInt)
    return hash(i(g), hash(j(g), hash(k(g), hash(α(g), hash(c(g), h)))))
end

function Base.:(==)(g_1::CartesianGaussian, g_2::CartesianGaussian)
    return i(g_1) == i(g_2) &&
           j(g_1) == j(g_2) &&
           k(g_1) == k(g_2) &&
           α(g_1) == α(g_2) &&
           c(g_1) == c(g_2)
end

"""
A Hermite Gaussian is

  ∂ⁱ/∂c₁ⁱ ∂ʲ/∂c₂ʲ ∂ᵏ/∂c₃ᵏ exp(-α(r - c)²)
= Λᵢ(x, α) Λⱼ(y, α) Λₖ(z, α) exp(-α(r - c)²)
"""
struct HermiteGaussian <: AbstractGaussian
    i::Int
    j::Int
    k::Int
    α::Float64
    c::Vector{Float64}
end

"""
A contracted Gaussian is just a sum of Gaussians.

cg(i, j, k) = ∑(ₙ, cₙ / |g(i, j, k, αₙ))| * g(i, j, k, αₙ))

The stored coefficients are normalized: cₙ / |g(i, j, k, αₙ))|

The number from BSE is cₙ, which is unnormalized.
"""
struct ContractedGaussian{T<:AbstractGaussian}
    coefficients::Vector{<:Number}
    gaussians::Vector{T}
end

"""
Normalized coefficients: cₙ / |g(i, j, k, αₙ))|
"""
coefficients(g::ContractedGaussian) = g.coefficients

function unnormalized_coefficients(g::ContractedGaussian)
    return coefficients(g) ./ normalization(gaussians(g))
end

"""
Unnormalized Gaussians: g(i, j, k, αₙ)
"""
gaussians(g::ContractedGaussian) = g.gaussians

Base.iterate(cg::ContractedGaussian) = ((first(coefficients(cg)), first(gaussians(cg))), 1)
function Base.iterate(cg::ContractedGaussian, i::Int)
    i >= length(coefficients(cg)) && return nothing
    ((coefficients(cg)[i+1], gaussians(cg)[i+1]), i + 1)
end

Base.length(cg::ContractedGaussian) = length(coefficients(cg))

"""
The contracted basis sets don't seem to come normalized out of the box. So one may need to normalize them.
"""
function contract(coefficients::Vector{<:Number}, gaussians::Vector{T}, normalize=false) where {T<:AbstractGaussian}
    cg = ContractedGaussian{T}(coefficients, gaussians)
    if normalize
        normalization = overlap_integral(gaussian_product(cg, cg))
        return sqrt(1 / normalization) * cg
    else
        return cg
    end
end

function combine(ccg::ContractedGaussian{CartesianGaussian})
    coeffs = coefficients(ccg)
    gs = gaussians(ccg)
    d = Dict{CartesianGaussian,Number}(g => 0.0 for g in gs)
    for (c, g) in zip(coeffs, gs)
        d[g] += c
    end

    filter!(((_, c),) -> abs(c) > 1e-10, d)

    return CCG(collect(values(d)), collect(keys(d)))
end

function contract(coefficients::Vector{Float64}, gaussians::Vector{HermiteGaussian})
    ContractedGaussian{HermiteGaussian}(coefficients, gaussians)
end

function Base.:+(cg_1::ContractedGaussian{T}, cg_2::ContractedGaussian{T}) where {T<:AbstractGaussian}
    contract(
        vcat(coefficients(cg_1), coefficients(cg_2)),
        vcat(gaussians(cg_1), gaussians(cg_2)),
    )
end

function Base.:*(c::Number, cg::ContractedGaussian{T}) where {T<:AbstractGaussian}
    contract(c * coefficients(cg), gaussians(cg))
end


i(g::AbstractGaussian)::Int = g.i
set_i(g::AbstractGaussian, new_i::Int) = Accessors.@set g.i = new_i
j(g::AbstractGaussian)::Int = g.j
set_j(g::AbstractGaussian, new_j::Int) = Accessors.@set g.j = new_j
k(g::AbstractGaussian)::Int = g.k
set_k(g::AbstractGaussian, new_k::Int) = Accessors.@set g.k = new_k
α(g::AbstractGaussian)::Float64 = g.α
set_α(g::AbstractGaussian, new_α::Int) = Accessors.@set g.α = new_α
c(g::AbstractGaussian)::Vector{Float64} = g.c
set_c(g::AbstractGaussian, new_c::Vector{Float64}) = Accessors.@set g.c = new_c


"""
Acronym for contracted Cartesian Gaussians.
"""
const CCG = ContractedGaussian{CartesianGaussian}

const CHG = ContractedGaussian{HermiteGaussian}
const Gaussian = Union{AbstractGaussian,ContractedGaussian}

struct Conj{T<:Gaussian}
    g::T
end

Base.adjoint(g::ContractedGaussian) = Conj(g)
Base.adjoint(g::T) where {T<:AbstractGaussian} = Conj{T}(g)
term(c::Conj) = c.g
