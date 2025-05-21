export spherical_to_cartesian, to_contracted_gaussian, SphericalGaussian

using SphericalHarmonics, Memoize

struct SphericalGaussian <: AbstractGaussian
    l::Int
    m::Int
    n::Int
    α::Float64
    c::Vector{Float64}
    gaussian_cache_number::Int
end

SphericalGaussian(l, m, n, α, c) = SphericalGaussian(l, m, n, α, c, -1)

get_l(g::SphericalGaussian) = g.l
get_m(g::SphericalGaussian) = g.m
get_n(g::SphericalGaussian) = g.n


function evaluate(g::SphericalGaussian, r::AbstractVector)
    x, y, z = r
    θ = acos(z / abs(r))
    ϕ = atan(y / x)

    P = computeYlm(θ, ϕ; lmax=get_l(g))
    return P[(get_l(g), get_m(g))] * exp(-α(g) * norm(r - c(g))^2)
end

"""
NOT Verified. Use at your own peril.

This does not have to work correctly for the integrals to be correct.
"""
function normalization(s::SphericalGaussian)
    n = get_l(s)
    #= return ((factorial(2n + 2) * π^(1 / 2)) /
            (2^(2n + 3) * factorial(n + 1) * α(s)^(n + 3 / 2)))^(-1 / 2) =#
    return 2^(-5 / 2 - n) * α(s)^(-3 / 2 - n) * gamma(3 / 2 + n)
end

function spherical_to_cartesian(l, m, l_x, l_y, l_z)
    m < 0 && return conj(spherical_to_cartesian(l, -m, l_x, l_y, l_z))
    t_1 = sqrt(factorial(2l_x) * factorial(2l_y) * factorial(2l_z) * factorial(l) * factorial(l - abs(m)) /
               (factorial(2l) * factorial(l_x) * factorial(l_y) * factorial(l_z) * factorial(l + abs(m)))) *
          1 / (2^l * factorial(l))
    j, r = divrem(l_x + l_y - abs(m), 2)
    r == 0 && j >= 0 || return 0

    t_2 = sum(i -> binomial(l, i) * binomial(i, j) * ((-1)^i * factorial(2l - 2i)) / factorial(l - abs(m) - 2i), 0:div(l - abs(m), 2))
    t_3 = sum(k -> binomial(j, k) * binomial(abs(m), l_x - 2k) * (complex(-1))^(sign(m) * (abs(m) - l_x + 2k) / 2), 0:j)
    return t_1 * t_2 * t_3
end

"""
Expand a unnormalized spherical Gaussians as a sum of unnormalized Cartesian Gaussians.
S(l,m) = ∑((r,s,t), |S|/|C| c(l,m,r,s,t) C(r,s,t))
"""
function to_contracted_gaussian(s::SphericalGaussian)
    l, m = get_l(s), get_m(s)
    lvs = filter(lv -> sum(lv) == l, vec(collect(Iterators.product(0:l, 0:l, 0:l))))

    coeffs = map(lv -> spherical_to_cartesian(l, m, lv...), lvs)
    gaussians = map(lv -> CartesianGaussian(lv..., α(s), c(s)), lvs)

return 1/normalization(s) * contract(coeffs .* normalization.(gaussians), gaussians)
end

"""
Expand a sum of unnormalized spherical Gaussians.
"""
function to_contracted_gaussian(cs::ContractedGaussian{SphericalGaussian})
    return combine(sum(coefficients(cs) .* to_contracted_gaussian.(gaussians(cs))))
end

to_contracted_gaussian(cs::AbstractGaussian) = cs
to_contracted_gaussian(cs::CCG) = cs

