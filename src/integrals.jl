using SpecialFunctions
using LinearAlgebra
using IterTools: product
using FastGaussQuadrature
using HypergeometricFunctions


export i,
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
    l,
    m,
    VNuc,
    kinetic_integral

"""
Compute the coefficients of the Hermite Gaussians for a product of two Cartesian Gaussians.
"""
function product_coefficients(n_a, n_b, N, α, pa, pb)::Float64
    d(n_a, n_b, N) = product_coefficients(n_a, n_b, N, α, pa, pb)
    0 <= N <= n_a + n_b || return 0
    n_a != 0 && return let n_a = n_a - 1
        1 / (2α) * d(n_a, n_b, N - 1) + pa * d(n_a, n_b, N) + (N + 1) * d(n_a, n_b, N + 1)
    end
    n_b != 0 && return let n_b = n_b - 1
        1 / (2α) * d(n_a, n_b, N - 1) + pb * d(n_a, n_b, N) + (N + 1) * d(n_a, n_b, N + 1)
    end

    return 1.0
end

Base.:*(a::CartesianGaussian, b::CartesianGaussian) = gaussian_product(a, b)

"""
The product of two cartesian Gaussians. The result is a contracted Hermite Gaussian.
"""
function gaussian_product(a::CartesianGaussian, b::CartesianGaussian)
    αp = α(a) + α(b)
    p = (c(a) * α(a) + c(b) * α(b)) / αp
    pa, pb = p - c(a), p - c(b)
    e = exp(-α(a) * α(b) / αp * norm(c(a) - c(b))^2)
    Ni, Nj, Nk = i(a) + i(b), j(a) + j(b), k(a) + k(b)
    DI = [product_coefficients(i(a), i(b), N, αp, pa[1], pb[1]) for N = 0:Ni]
    DJ = [product_coefficients(j(a), j(b), N, αp, pa[2], pb[2]) for N = 0:Nj]
    DK = [product_coefficients(k(a), k(b), N, αp, pa[3], pb[3]) for N = 0:Nk]

    coefficients_and_gaussians = [
        ((e * DI[ni+1] * DJ[nj+1] * DK[nk+1]), HermiteGaussian(ni, nj, nk, αp, p)) for
        ni = 0:Ni, nj = 0:Nj, nk = 0:Nk
    ] |> vec

    contract(first.(coefficients_and_gaussians), last.(coefficients_and_gaussians))
end

Base.:*(a::ContractedGaussian{CartesianGaussian}, b::ContractedGaussian{CartesianGaussian}) = gaussian_product(a, b)

"""
The product of two contracted Cartesian Gaussians.
"""
function gaussian_product(a::ContractedGaussian{CartesianGaussian},
    b::ContractedGaussian{CartesianGaussian})
    contracted_gaussians = [
        (c_a * c_b * gaussian_product(g_a, g_b))
        for (c_a, g_a) in zip(coefficients(a), gaussians(a)),
        (c_b, g_b) in zip(coefficients(b), gaussians(b))
    ] |> vec
    sum(contracted_gaussians)
end

# Overlap integral of a Hermite Gaussian, which arises from the product of two Gaussians.


"""
The overlap integral of two Hermite Gaussians
"""
function overlap_integral(g::HermiteGaussian)
    (i(g) == 0 && j(g) == 0 && k(g) == 0) || return 0
    return (π / α(g))^(3 / 2)
end

Base.:*(a::Conj{T}, b::T) where {T} = overlap_integral(gaussian_product(term(a), b))

"""
The overlap integral of two contracted Hermite Gaussians.
"""
function overlap_integral(cg::ContractedGaussian{HermiteGaussian})
    # sum(zip(coefficients(cg), gaussians(cg))) do (c, g)
    #     c * overlap_integral(g)
    # end
    sum(cg) do (c, g)
        c * overlap_integral(g)
    end

end


function kinetic_integral(p::ContractedGaussian, q::ContractedGaussian)
    sum(zip(coefficients(p), gaussians(p))) do (c_p, g_p)
        sum(zip(coefficients(q), gaussians(q))) do (c_q, g_q)
            c_p * c_q * kinetic_integral(g_p, g_q)
        end
    end
end

function kinetic_integral(p::CartesianGaussian, q::CartesianGaussian)
    s = p' * q
    sum(zip([i, j, k], [set_i, set_j, set_k])) do (x, set_x)
        α(q) * (2x(q) + 1) * s - 2α(q)^2 * (p' * set_x(q, x(q) + 2)) -
        0.5x(q) * (x(q) - 1) * (x(q) >= 2 && (p' * set_x(q, x(q) - 2)))
    end
end


function two_electron_integral(p::ContractedGaussian{HermiteGaussian}, q::ContractedGaussian{HermiteGaussian})
    sum(zip(coefficients(p), gaussians(p))) do (c_p, g_p)
        sum(zip(coefficients(q), gaussians(q))) do (c_q, g_q)
            c_p * c_q * two_electron_integral(g_p, g_q)
        end
    end
end

function two_electron_integral(p::HermiteGaussian, q::HermiteGaussian)
    λ = 2 * π^(5 / 2) * α(p)^(-1) * α(q)^(-1) * (α(p) + α(q))^(-1 / 2)
    αc = α(p) * α(q) / (α(p) + α(q))
    pq = c(p) - c(q)
    λ * (-1)^(i(q) + j(q) + k(q)) * R_tensor(i(p) + i(q), j(p) + j(q), k(p) + k(q), 0, αc, pq...)
end


function evaluate(g::CartesianGaussian, r::Vector{<:Real})
    x, y, z = r
    c1, c2, c3 = c(g)
    (x - c1)^i(g) * (y - c2)^j(g) * (z - c3)^k(g) * exp(-α(g) * norm(r - c(g))^2)
end

function evaluate(g::HermiteGaussian, r::Vector{<:Real})
    x, y, z = r
    c1, c2, c3 = c(g)
    Lambda(i(g), x - c1, α(g)) * Lambda(j(g), y - c2, α(g)) * Lambda(k(g), z - c3, α(g)) * exp(-α(g) * norm(r - c(g))^2)
end

function evaluate(g::ContractedGaussian, r::Vector{<:Real})
    sum(zip(coefficients(g), gaussians(g))) do (c, g)
        c * evaluate(g, r)
    end
end

function hermite(N::Int, x::Number)
    N == 0 && return 1
    N == 1 && return 2x
    return 2x * hermite(N - 1, x) - 2 * (N - 1) * hermite(N - 2, x)
end

function Lambda(N::Int, x, α)
    α^(N / 2) * hermite(N, α^(1 / 2) * x)
end

function double_factorial(n::Int)
    n <= 1 && return 1
    n * double_factorial(n - 2)
end


"""
The normalization factor of a Cartesian Gaussian.

lookup the result from cache instead if a cache is given.
"""
function normalization(g::CartesianGaussian)
    # haskey(cache.cartesian_normalization_factors, cache_number(g)) &&
    #     return cache.cartesian_normalization_factors[cache_number(g)]
    n(v) = (2α(g) / π)^(1 / 4) * (4α(g))^(v / 2) * (double_factorial(2v - 1))^(-1 / 2)
    prod(n, [i(g), j(g), k(g)])
end

# function normalization(g::SphericalGaussian)
#     n = l(g)
#     (factorial(2n + 2) * π^(1/2) /
#      (2^(2n+3) * factorial(n+1) * α^(n+3/2)))^(-1/2)
# end

struct VNuc
    atoms::Vector{<:AbstractAtom}
end

Base.:+(v_1::VNuc, v_2::VNuc) = VNuc(vcat(atoms(v_1), atoms(v_2)))

atoms(v::VNuc) = v.atoms

struct VNucState{T<:Union{Gaussian,Conj}}
    v_nuc::VNuc
    state::T
end

atoms(v::VNucState) = atoms(v.v_nuc)
state(v::VNucState) = v.state


Base.:*(v::VNuc, g::T) where {T<:Gaussian} = VNucState{T}(v, g)
Base.:*(g::T, v::VNuc) where {T<:Conj} = VNucState{T}(v, g)
Base.:*(v::VNucState{T}, g::S) where {T<:Conj,S<:Gaussian} =
    sum(atoms(v)) do a
        nuclear_potential(term(state(v)) * g, a)
    end
Base.:*(g::S, v::VNucState{T}) where {S<:Conj,T<:Gaussian} =
    sum(atoms(v)) do a
        nuclear_potential(term(g) * state(v), a)
    end


function nuclear_potential(cg::ContractedGaussian{HermiteGaussian}, a::AbstractAtom)
    -charge(a) * sum(zip(coefficients(cg), gaussians(cg))) do (c, g)
        c * nuclear_potential(g, coordinates(a))
    end
end

# function nuclear_potential(cg_1::ContractedGaussian{CartesianGaussian}, cg_2::ContractedGaussian{CartesianGaussian}, a::AbstractAtom)
#     v = 0
#     for (c_1, g_1) in zip(coefficients(cg_1), gaussians(cg_1))
#         for (c_2, g_2) in zip(coefficients(cg_2), gaussians(cg_2))
#             v += c_1 * c_2 * nuclear_potential(g_1 * g_2, a)
#         end
#     end
#     return v
# end

"""
The sign of `c(g) -r` is not explained in MD and I didn't pay attention to it.
It turn into a bug that took me a day to fix.
"""
function nuclear_potential(g::HermiteGaussian, r::Vector{<:Real})
    return (2π / α(g)) * R_tensor(i(g), j(g), k(g), 0, α(g), (c(g) - r)...)
end

function R_tensor(i::Int, j::Int, k::Int, n::Int, α::Float64, a::Float64, b::Float64, c::Float64)
    (i < 0 || j < 0 || k < 0) && return 0
    T = α * (a^2 + b^2 + c^2)
    (i == j == k == 0) && return (-2α)^n * F(n, T)
    R(i::Int, j::Int, k::Int, n::Int) = R_tensor(i, j, k, n, α, a, b, c)
    i == j == 0 && return c * R(0, 0, k - 1, n + 1) + (k - 1) * R(0, 0, k - 2, n + 1)
    i == 0 && return b * R(0, j - 1, k, n + 1) + (j - 1) * R(0, j - 2, k, n + 1)
    return a * R(i - 1, j, k, n + 1) + (i - 1) * R(i - 2, j, k, n + 1)
end


"""
This is currently not properly done.
The function is essential Γ(1/2 + n, T) / T^(1/2+n).
The problem with this is the fake singularity at T = 0.
THe current workaround is to Taylor expand near 0.

Other ways of evaluating this includes the Hyper Geometric functions,
but the one provided in HypergeometricFunctions.jl has not worked well here.

One can also just do a quadrature, but hundreds of quadrature points are needed
to converge.
"""
function F(n::Int, T::Float64)
    # x, w = gausslegendre(100 + 2n)
    # f(u) = u^(2n) * exp(-T * u^2)
    # return dot(w, f.(x)) / 2
    # HypergeometricFunctions.M(n+1.5, n+0.5, -T) / (2n + 1)

    if abs(T) < 0.01
        return 1 / (1 + 2n) +
               T / (-3 - 2n) +
               T^2 / (2 * (5 + 2n)) -
               T^3 / (12 * (7 / 2 + n)) +
               T^4 / (48 * (9 / 2 + n)) -
               T^5 / (240 * (11 / 2 + n)) +
               T^6 / (1440 * (13 / 2 + n)) -
               T^7 / (10080 * (15 / 2 + n))
    else
        return 1 / 2 * gamma(1 / 2 + n) * first(gamma_inc(1 / 2 + n, T)) / T^(1 / 2 + n)
    end
end


    # t_0 = α(q) * (2 * (i(q) + j(q) + k(q)) + 3) * (p' * q)
    # t_1 = -2 * α(q)^2 * (p' * set_i(q, i(q)+2) +
    #                      p' * set_j(q, j(q)+2) +
    #                      p' * set_k(q, k(q)+2))
    # t_2 = -0.5 * (i(q) * (i(q) - 1) * (i(q) >= 2 && (p' * set_i(q, i(q) - 2))) + 
    #               j(q) * (j(q) - 1) * (j(q) >= 2 && (p' * set_j(q, j(q) - 2))) +
    #               k(q) * (k(q) - 1) * (k(q) >= 2 && (p' * set_k(q, k(q) - 2))))

    # return t_0 + t_1 + t_2