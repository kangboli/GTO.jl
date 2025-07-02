using SpecialFunctions
using LinearAlgebra
using HypergeometricFunctions
# using ApproxFun


export
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
    VNuc,
    kinetic_integral,
    ∇²,
    pos_integral,
    r2_integral,
    normalization,
    @r_unroll

"""
Compute the coefficients of the Hermite Gaussians for a product of two Cartesian Gaussians.
"""
function product_coefficients(n_a::Int, n_b::Int, N::Int, α::Float64, pa::Float64, pb::Float64)::Float64
    d(n_a::Int, n_b::Int, N::Int) = product_coefficients(n_a, n_b, N, α, pa, pb)
    0 <= N <= n_a + n_b || return 0
    n_a != 0 && return let n_a = n_a - 1
        d(n_a, n_b, N - 1) / (2α) + pa * d(n_a, n_b, N) + (N + 1) * d(n_a, n_b, N + 1)
    end
    n_b != 0 && return let n_b = n_b - 1
        d(n_a, n_b, N - 1) / (2α) + pb * d(n_a, n_b, N) + (N + 1) * d(n_a, n_b, N + 1)
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
    e = exp(-α(a) * α(b) * norm(c(a) - c(b))^2 / αp)
    Ni, Nj, Nk = i(a) + i(b), j(a) + j(b), k(a) + k(b)
    DI = [product_coefficients(i(a), i(b), N, αp, pa[1], pb[1]) for N = 0:Ni]
    DJ = [product_coefficients(j(a), j(b), N, αp, pa[2], pb[2]) for N = 0:Nj]
    DK = [product_coefficients(k(a), k(b), N, αp, pa[3], pb[3]) for N = 0:Nk]
    coefficients = Vector{Float64}()
    hermites = Vector{HermiteGaussian}()
    for ni in 0:Ni, nj in 0:Nj, nk in 0:Nk
        coeff = e * DI[ni+1] * DJ[nj+1] * DK[nk+1]
        #= abs(coeff) < 1e-12 && continue =#
        push!(hermites, HermiteGaussian(ni, nj, nk, αp, p))
        push!(coefficients, coeff)
    end

    contract(coefficients, hermites)
end

Base.:*(a::CCG, b::CCG) = gaussian_product(a, b)

"""
The product of two contracted Cartesian Gaussians.
"""
function gaussian_product(a::CCG, b::CCG)
    [(c_a' * c_b * gaussian_product(g_a, g_b)) for (c_a, g_a) in a, (c_b, g_b) in b] |> vec |> sum
end

# Overlap integral of a Hermite Gaussian, which arises from the product of two Gaussians.


"""
The overlap integral of a single Hermite Gaussian.
"""
function overlap_integral(g::HermiteGaussian)
    (i(g) == 0 && j(g) == 0 && k(g) == 0) || return 0
    return (π / α(g))^(3 / 2)
end

Base.:*(a::Conj{T}, b::T) where {T} = overlap_integral(term(a) * b)

"""
The overlap integral of two Cartesian Gaussians.
"""
overlap_integral(cg::CHG) = sum(((c, g),) -> c * overlap_integral(g), cg)

struct Kinetic end

const ∇² = Kinetic()
struct KineticState{T<:Union{Gaussian,Conj}}
    state::T
end

state(k::KineticState) = k.state

Base.:*(::Kinetic, g::T) where {T<:Gaussian} = KineticState{T}(g)
Base.:*(g::T, ::Kinetic) where {T<:Conj} = KineticState{T}(g)
Base.:*(k::KineticState{T}, g::S) where {T<:Conj,S<:Gaussian} = kinetic_integral(term(state(k)), g)
Base.:*(g::S, k::KineticState{T}) where {T<:Gaussian,S<:Conj} = kinetic_integral(term(g), state(k))



function kinetic_integral(p::ContractedGaussian, q::ContractedGaussian)
    sum(p) do (c_p, g_p)
        sum(((c_q, g_q),) -> c_p * c_q * kinetic_integral(g_p, g_q), q)
    end
end

function kinetic_integral(p::CartesianGaussian, q::CartesianGaussian)
    s = p' * q
    sum(zip([i, j, k], [set_i, set_j, set_k])) do (x, set_x)
        α(q) * (2x(q) + 1) * s - 2α(q)^2 * (p' * set_x(q, x(q) + 2)) -
        0.5x(q) * (x(q) - 1) * (x(q) >= 2 && (p' * set_x(q, x(q) - 2)))
    end
end

Base.:|(p::CHG, q::CHG) = two_electron_integral(p, q)

function two_electron_integral(p::CHG, q::CHG)
    sum(p) do (c_p, g_p)
        sum(((c_q, g_q),) -> c_p' * c_q * two_electron_integral(g_p, g_q), q)
    end
end

function two_electron_integral(p::HermiteGaussian, q::HermiteGaussian)
    αp, αq = α(p), α(q)
    λ = 2 * π^(5 / 2) * αp^(-1) * αq^(-1) * (αp + αq)^(-1 / 2)
    αc = αp * αq / (αp + αq)
    pq = p.c - q.c
    λ * (-1)^(i(q) + j(q) + k(q)) * R_tensor(i(p) + i(q), j(p) + j(q), k(p) + k(q), 0, αc, pq...)
end


function evaluate(g::CartesianGaussian, r::AbstractVector)
    x, y, z = r
    c1, c2, c3 = c(g)
    (x - c1)^i(g) * (y - c2)^j(g) * (z - c3)^k(g) * exp(-α(g) * norm(r - c(g))^2)
end

function evaluate(g::HermiteGaussian, r::AbstractVector)
    x, y, z = r
    c1, c2, c3 = c(g)
    Lambda(i(g), x - c1, α(g)) * Lambda(j(g), y - c2, α(g)) *
    Lambda(k(g), z - c3, α(g)) * exp(-α(g) * norm(r - c(g))^2)
end

evaluate(g::ContractedGaussian, r::AbstractVector) = sum(((c, g),) -> c * evaluate(g, r), g)

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
The normalization factor (1/|g|) of a Cartesian Gaussian.

Note that this is the inverse of the norm.
g̃ = 1/|g| g

"""
function normalization(g::CartesianGaussian)
    # haskey(cache.cartesian_normalization_factors, cache_number(g)) &&
    #     return cache.cartesian_normalization_factors[cache_number(g)]
    n(v) = (2α(g) / π)^(1 / 4) * (4α(g))^(v / 2) * (double_factorial(2v - 1))^(-1 / 2)
    prod(n, [i(g), j(g), k(g)])
end

function normalization(cg::CHG)
    sum(cg) do (c_1, g_1)
        sum(cg) do (c_2, g_2)
            i(g_1) == j(g_1) == k(g_1) == i(g_2) == j(g_2) == k(g_2) == 0 ?
            c_1 * c_2 * (π / (α(g_1) + α(g_2)))^(3 / 2) : 0
        end
    end

end

struct VNuc
    atoms::Vector{<:AbstractAtom}
end

VNuc(atoms::Vararg{<:AbstractAtom}) = VNuc(collect(atoms))

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
    sum(a -> nuclear_potential(term(state(v)) * g, a), atoms(v))

Base.:*(g::S, v::VNucState{T}) where {S<:Conj,T<:Gaussian} =
    sum(a -> nuclear_potential(term(g) * state(v), a), atoms(v))


function nuclear_potential(cg::CHG, a::AbstractAtom)
    -charge(a) * sum(((c, g),) -> c * nuclear_potential(g, coordinates(a)), cg)
end

"""
The sign of `c(g) -r` is not explained in MD and I didn't pay attention to it.
It turn into a bug that took me a day to fix.
"""
function nuclear_potential(g::HermiteGaussian, r::AbstractVector{<:Real})::Float64
    return (2π / α(g)) * R_tensor(i(g), j(g), k(g), 0, α(g), (c(g) - r)...)
end

function R_tensor(i::Int, j::Int, k::Int, n::Int, α::Float64, a::Float64, b::Float64, c::Float64, T::Float64=α * (a^2 + b^2 + c^2))::Float64
    (i < 0 || j < 0 || k < 0) && return 0.0
    (i == j == k == 0) && return (-2α)^n * F(n, T)
    i == j == 0 && return c * R_tensor(0, 0, k - 1, n + 1, α, a, b, c, T) + (k - 1) * R_tensor(0, 0, k - 2, n + 1, α::Float64, a, b, c, T)
    i == 0 && return b * R_tensor(0, j - 1, k, n + 1, α, a, b, c, T) + (j - 1) * R_tensor(0, j - 2, k, n + 1, α::Float64, a, b, c, T)
    return a * R_tensor(i - 1, j, k, n + 1, α, a, b, c, T) + (i - 1) * R_tensor(i - 2, j, k, n + 1, α::Float64, a, b, c, T)
end

l_max = 14
generated_integrals = Array{Function,3}(undef, l_max, l_max, l_max)

#= function R_tensor(i::Int, j::Int, k::Int, n::Int, α::Float64, a::Float64, b::Float64, c::Float64, T::Float64=α * (a^2 + b^2 + c^2))::Float64
    @assert n == 0
    GTO.generated_integrals[i+1, j+1, k+1](α, a, b, c, T)
    #= GTO.gen_int(Val(i+1), Val(j+1), Val(k+1), α, a, b, c, T) =#
end =#

macro r_unroll(i, j, k, n, α, a, b, c, T)
    (i < 0 || j < 0 || k < 0) && return esc(:(0.0))
    (i == j == k == 0) && return esc(:((-2 * $α)^$n * F($n, $T)))

    if i == j == 0
        k_term = :($c * @r_unroll(0, 0, $(k - 1), $(n + 1), $α, $a, $b, $c, $T))
        k > 1 || return esc(k_term)
        return esc(:($k_term + $(k - 1) * @r_unroll(0, 0, $(k - 2), $(n + 1), $α, $a, $b, $c, $T)))
    end

    if i == 0
        j_term = :(b * @r_unroll(0, $(j - 1), $k, $(n + 1), $α, $a, $b, $c, $T))
        j > 1 || return esc(j_term)
        return esc(:($j_term + $(j - 1) * @r_unroll(0, $(j - 2), $k, $(n + 1), $α, $a, $b, $c, $T)))
    end

    i_term = :($a * @r_unroll($(i - 1), $j, $k, $(n + 1), $α, $a, $b, $c, $T))

    i > 1 || return esc(i_term)
    return esc(:($i_term + $(i - 1) * @r_unroll($(i - 2), $j, $k, $(n + 1), $α, $a, $b, $c, $T)))
end

for i in 1:l_max, j in 1:l_max, k in 1:l_max
    i + j + k <= l_max || continue
    generated_integrals[i, j, k] = eval(:((α::Float64, a::Float64, b::Float64, c::Float64, T::Float64) -> @r_unroll($(i - 1), $(j - 1), $(k - 1), 0, α, a, b, c, T)))
end

function gen_int()
    0.0
end

for i in 1:l_max, j in 1:l_max, k in 1:l_max
    i + j + k <= l_max || continue
    eval(:(
        function GTO.gen_int(::Val{$i}, ::Val{$j}, ::Val{$k}, α::Float64, a::Float64, b::Float64, c::Float64, T::Float64)
            @r_unroll($(i - 1), $(j - 1), $(k - 1), 0, α, a, b, c, T)
        end
    ))
end


"""
The Boys function
"""
function F(n::Int, T::Number)::Float64
    #= abs(T) < 0.01 && return 1 / (1 + 2n) +
                            T / (-3 - 2n) +
                            T^2 / (2 * (5 + 2n)) -
                            T^3 / (12 * (7 / 2 + n)) +
                            T^4 / (48 * (9 / 2 + n)) -
                            T^5 / (240 * (11 / 2 + n)) +
                            T^6 / (1440 * (13 / 2 + n)) -
                            T^7 / (10080 * (15 / 2 + n))

    m = 1 / 2 + n
    return 1 / 2 * (gamma(m) * first(gamma_inc(m, T))) / T^(m) =#
    _₁F₁(n + 0.5, n + 1.5, -T) / (2.0 * n + 1.0)
end

function pos_integral(g_1::CartesianGaussian, g_2::CartesianGaussian)
    s = overlap_integral(g_1 * g_2)
    r = c(g_2) * s
    g_i = set_i(g_2, i(g_2) + 1)
    g_j = set_j(g_2, j(g_2) + 1)
    g_k = set_k(g_2, k(g_2) + 1)
    r[1] += overlap_integral(g_1 * g_i)
    r[2] += overlap_integral(g_1 * g_j)
    r[3] += overlap_integral(g_1 * g_k)
    return r
end

function pos_integral(g_1::CCG, g_2::CCG)
    sum(g_1) do (c_1, g_1)
        sum(g_2) do (c_2, g_2)
            c_1' * c_2 * pos_integral(g_1, g_2)
        end
    end
end

function r2_integral(g_1::CartesianGaussian, g_2::CartesianGaussian)
    x = c(g_2)
    integral = -norm(x)^2 * overlap_integral(g_1 * g_2)
    integral += 2dot(x, pos_integral(g_1, g_2))
    g_ii = set_i(g_2, i(g_2) + 2)
    g_jj = set_j(g_2, j(g_2) + 2)
    g_kk = set_k(g_2, k(g_2) + 2)
    r2_g = sum(g -> overlap_integral(g_1 * g), [g_ii, g_jj, g_kk])
    integral += r2_g
    return integral
end

function r2_integral(g_1::CCG, g_2::CCG)
    sum(g_1) do (c_1, g_1)
        sum(g_2) do (c_2, g_2)
            c_1' * c_2 * r2_integral(g_1, g_2)
        end
    end
end



#= struct Pos end

struct PosState{T<:Union{Gaussian,Conj}}
    state::T
end

Base.:*(::Pos, state::T) where T <: Gaussian = PosState{T}(state)
Base.:*(state::T, ::Pos) where T <: Conj = PosState{T}(state)

function Base.:*(v::PosState{T}, g::S) where {T<:Conj,S<:Gaussian} 
end

function Base.:*(g::S, v::PosState{T}) where {S<:Conj,T<:Gaussian} 
end

 =#
