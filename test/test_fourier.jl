using GTO, WTP, Test

b = load_basis("sto-3g.0.json")

h = make_atom(:H, 1, 0, 0)

basis = make_gaussians(b, h)
g = basis[1][0][1].gaussians[1]

ĝ = ft_gaussian(g)
ĝ_0 = ft_gaussian_simplified(g)

@test isapprox(ĝ([1, 1, 1]), ĝ_0([1, 1, 1]), atol=1e-10)


h = make_atom(:C, 0, 0, 0)

basis = make_gaussians(b, h)
g = basis[2][1][1].gaussians[1]

ĝ = ft_gaussian(g)
ĝ_0 = ft_gaussian_simplified(g)

@test isapprox(ĝ([1, 0, 0]), ĝ_0([1, 0, 0]), atol=1e-10)


ĝ_num = ft_numerical(g)
@test isapprox(ĝ([1, 0, 0]), ĝ_num([1, 0, 0]), atol=1e-2)

d = hermite_coefficient(3)
@test d[3] == 8
@test d[1] == -12

d = hermite_coefficient(4)
@test d[4] == 16
@test d[2] == -48
@test d[0] == 12

hermite_factors(4, 0.3)

h = hermite_to_cart(gaussians(g * g)[3])
@test isapprox(h.coefficients, [-11.765, 138.415], atol=1e-3)



