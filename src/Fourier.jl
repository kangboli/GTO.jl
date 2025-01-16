export ft_single, ft_single_0, ft_gaussian, ft_gaussian_simplified, ft_numerical, hermite_coefficient, hermite_factors, hermite_to_cart

"""
Fourier transform of Gaussian type orbitals.
"""

function ft_single(l, c, α)
    return G -> begin
        1 / 2 * exp(-α * c^2)α^(-l / 2) * (
            (1 + (-1)^l) * gamma((1 + l) / 2) *
            _₁F₁((l + 1) / 2, 1 / 2, -(G + 2im * c * α)^2 / (4α)) / √(α) +
            im * (-1 + (-1)^l) * (G + 2im * c * α) * gamma(1 + l / 2) *
            _₁F₁((2 + l) / 2, 3 / 2, -(G + 2im * c * α)^2 / (4α)) / α)
    end
end

function ft_single_0(c, α)
    return G -> √(π) * exp(-α * c^2 - (G + 2im * α * c)^2 / (4α)) / √(α)
end

function ft_single_1(c, α)
    return G -> -(im * exp(-α * c^2 - (G + 2im * c * α)^2 / (4α)) * √(π) * (G + 2im * c * α)) / (2 * α^(3 / 2))
end

function ft_single_2(c, α)
    return G -> (exp(-α * c^2 - (G + 2im * c * α)^2 / (4α)) * √(π) * (1 - (G + 2im * c * α)^2 / (2α))) / (2 * α^(3 / 2))
end

function ft_gaussian(g::CartesianGaussian)
    return G -> begin
        ft_single(i(g), c(g)[1], α(g))(G[1]) *
        ft_single(j(g), c(g)[2], α(g))(G[2]) *
        ft_single(k(g), c(g)[3], α(g))(G[3])
    end
end

ft_simplified = [ft_single_0, ft_single_1, ft_single_2]

function ft_gaussian_simplified(g::CartesianGaussian)
    return G -> begin
        ft_simplified[i(g)+1](c(g)[1], α(g))(G[1]) *
        ft_simplified[j(g)+1](c(g)[2], α(g))(G[2]) *
        ft_simplified[k(g)+1](c(g)[3], α(g))(G[3])
    end
end


function ft_numerical(g::CartesianGaussian)
    domain = (-5, 4)
    n_grid = domain[2] - domain[1] + 1
    homecell = make_grid(HomeCell3D, tuple([6 * b / n_grid for b in CARTESIAN_BASIS]...), (domain, domain, domain))
    gaussian_on_grid = map(r -> evaluate(g, coordinates(r)), homecell)
    return G -> begin
        exp_on_grid = map(r -> exp(im * G' * coordinates(r)), homecell)
        return vec(elements(exp_on_grid))' * vec(elements(gaussian_on_grid)) * det(basis_matrix(homecell))
    end
end

function hermite_coefficient(n)
    iseven(n) && return Dict(
        2l => factorial(n) * (-1)^(n / 2 - l) / (
            factorial(2l) * factorial(div(n, 2) - l)) * 2^(2l)
        for l in 0:div(n, 2))
    return Dict(
        2 * l + 1 => factorial(n) * (-1)^((n - 1) / 2 - l) / (
            factorial(2l + 1) * factorial(div(n - 1, 2) - l)) * 2^(2l + 1)
        for l in 0:div(n - 1, 2))
end

function hermite_factors(n::Int, α)
    d = hermite_coefficient(n)
    for k in keys(d)
        d[k] *= α^(n / 2) * √(α)^(k)
    end
    return d
end

function hermite_to_cart(hg::HermiteGaussian)
    d_i = hermite_factors(i(hg), α(hg))
    d_j = hermite_factors(j(hg), α(hg))
    d_k = hermite_factors(k(hg), α(hg))
    coeff = Vector{Float64}()
    carts = Vector{CartesianGaussian}()
    for (k_i, k_j, k_k) in product(collect(keys(d_i)), collect(keys(d_j)), collect(keys(d_k)))
        cart_coeff = d_i[k_i] * d_j[k_j] * d_k[k_k]

        push!(carts, CartesianGaussian(k_i, k_j, k_k, α(hg), c(hg)))
        push!(coeff, cart_coeff)
    end

    return CCG(coeff, carts)
end
