using GTO
using Test

using Makie
using GLMakie


c1 = CartesianGaussian(0,0,0, 1, [0,0,0])
c2 = CartesianGaussian(0,0,0, 1, [1,0,0])
gaussian_product(c1, c2)



@test evaluate(c1, [1,1,0]) * evaluate(c2, [1,1,0]) == evaluate(gaussian_product(c1, c2), [1,1,0])


c1 = CartesianGaussian(1,0,0, 1, [0,0,0])
c2 = CartesianGaussian(0,0,0, 1, [1,0,0])
gaussian_product(c1, c2)

@test evaluate(c1, [0.5,1,1]) * evaluate(c2, [0.5,1,1]) == evaluate(gaussian_product(c1, c2), [0.5,1,1])

lines([evaluate(c1, [x,1,1]) * evaluate(c2, [x,1,1]) for x = 0:0.1:1])
lines!([evaluate(gaussian_product(c1, c2), [x,1,1]) for x = 0:0.1:1])

evaluate(gaussian_product(c1, c2), [0,1,1])
evaluate(c1, [0,1,1])

gaussian_product(c1, c2)

@test hermite(3,3) == 8 * 3^3 - 12 * 3



c1 = CartesianGaussian(1,3,1, 1, [0,0,1])
c2 = CartesianGaussian(0,2,1, 1, [0,1,0])
c3= gaussian_product(c1, c2)
isapprox(evaluate(gaussian_product(c1, c2), [1.2,.2,.3]), evaluate(c1, [1.2,.2,.3]) * evaluate(c2, [1.2,.2,.3]), atol=1e-6)

nuclear_potential(gaussians(c3)[1] , [2,2,2])


s = SphericalGaussian(3,1, 1, [0,0,0])

spherical_to_cartesian(s)

b = load_basis("basis_set_bundle-json-bib/sto-3g.1.json")

c = Cache()

a = Atom(6, [0,0,0], c)

shells = make_gaussians(b, a)
(shells |> last)[1]
