using Test, GTO, WTP

dir = "basis_set_bundle-json-bib"
for item in readdir(dir)
    endswith(item, "json") || continue
    occursin("cren", item) && continue
    occursin("ecp", item) && continue
    b = load_basis(item)

    for e in values(b.elements)
        for s in e.electron_shells
            #= @assert !(s.function_type == "gto" && any(i->i>1, s.angular_momentum)) =#
            #= @assert !(s.function_type == "gto_cartesian" && any(i->i<2, s.angular_momentum)) =#
            if s.function_type == "gto_spherical" && any(i -> i < 2, s.angular_momentum)
                println(s)
            end
            #= @assert !(s.function_type == "gto_cartesian" && any(i->i<2, s.angular_momentum)) =#
        end
    end

end

b = load_basis("cc-pvqz.0.json")
b = load_basis("sto-3g.0.json")
b = load_basis("6-31g.0.json")

spherical_to_cartesian(0, 0, 0, 0, 0) ≈ 1

spherical_to_cartesian(1, 0, 0, 0, 1) ≈ 1
spherical_to_cartesian(1, 1, 1, 0, 0) ≈ 1 / sqrt(2)
spherical_to_cartesian(1, 1, 0, 1, 0) ≈ im / sqrt(2)

spherical_to_cartesian(2, 0, 0, 0, 2) ≈ 1
spherical_to_cartesian(2, 0, 2, 0, 0) ≈ -1 / 2
spherical_to_cartesian(2, 0, 0, 2, 0) ≈ -1 / 2
spherical_to_cartesian(2, 1, 1, 0, 1) ≈ 1 / sqrt(2)
spherical_to_cartesian(2, 1, 0, 1, 1) ≈ im / sqrt(2)

spherical_to_cartesian(2, 2, 2, 0, 0) ≈ sqrt(3 / 8)
spherical_to_cartesian(2, 2, 0, 2, 0) ≈ -sqrt(3 / 8)
spherical_to_cartesian(2, 2, 1, 1, 0) ≈ im / sqrt(2)

spherical_to_cartesian(3, 0, 0, 0, 3) ≈ 1
spherical_to_cartesian(3, 0, 2, 0, 1) ≈ -3 / (2sqrt(5))
spherical_to_cartesian(3, 0, 0, 2, 1) ≈ -3 / (2sqrt(5))
spherical_to_cartesian(3, 1, 1, 0, 2) ≈ sqrt(3 / 5)
spherical_to_cartesian(3, 1, 0, 1, 2) ≈ im * sqrt(3 / 5)
spherical_to_cartesian(3, 1, 3, 0, 0) ≈ -sqrt(3) / 4
spherical_to_cartesian(3, 1, 0, 3, 0) ≈ -im * sqrt(3) / 4

spherical_to_cartesian(3, 1, 1, 2, 0) ≈ -sqrt(3) / (4sqrt(5))
spherical_to_cartesian(3, 1, 2, 1, 0) ≈ -im * sqrt(3) / (4sqrt(5))

sto3g = load_basis("sto-3g.1.json")
sto3g

zn = make_atom(:Zn, 0, 0, 0)

gaussians = generate_basis(sto3g, zn)

cartesians = to_contracted_gaussian.(gaussians)


S = [p' * q for p in cartesians, q in cartesians]
S = [p' * VNuc(zn) * q for p in cartesians, q in cartesians]
J = [(p * q | r * s) for p in cartesians, q in cartesians, r in cartesians, s in cartesians]
@time J = [(cartesians[1] * cartesians[1] | r * s) for r in cartesians, s in cartesians]
map(i -> abs(i) > 1e-10, J)
evaluate(cartesians[end] * cartesians[end-1], [0.2, 0.1, 0.1])
evaluate(cartesians[end], [0.1, 0.1, 0.02])
d = cartesians[end]
(cartesians[1] * cartesians[1] | cartesians[1] * cartesians[1])
d2 = d * d
@time (d * d | d * d)
f = (α, a, b, c, T) -> @r_unroll(1, 1, 1, 0, α, a, b, c, T)
f = (α, a, b, c, T) -> (a * (b * (c * ((-2α)^3 * F(3, T)))))
@macroexpand @r_unroll(1, 1, 1, 0, α, a, b, c, T)
@time f(1, 1, 1, 1, 3)
@time R_tensor(1, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 3.0)

d2.coefficients[end]
(@elapsed (d2 | d2)) * 18^4
length(d2.gaussians)
t = @elapsed two_electron_integral(d2.gaussians[2], d2.gaussians[10])
t * 36^2

gaussians = make_gaussians(sto3g, zn)
basis = vcat(get_basis.(gaussians)...)
cartesians = to_contracted_gaussian.(vcat(get_basis.(gaussians)...))
basis[end]
d = to_contracted_gaussian(basis[end])
d' * d


S = [p' * q for p in cartesians, q in cartesians]
cartesians[end]' * cartesians[end]

basis[1]' * basis[1]

sg = SphericalGaussian(0, 0, 0, 1, [0, 0, 0])
ccg = to_contracted_gaussian(sg)
sqrt(ccg' * ccg) ≈ 1 / normalization(sg)


cg = CartesianGaussian(0, 0, 0, 1, [0, 0, 0])

normalization(sg) / normalization(cg)


0.2821
normalization(cg)
1 / π^(3 / 4)

normalization(cg)
normalization(sg) / normalization(cg)

2 / π^(1 / 4)
sqrt(cg' * cg) * normalization(cg)

1 / (2sqrt(π)) / normalization(sg)^2
(1 / (2sqrt(π))) * (2 / (π^(1 / 4)))^2

normalization(cg) * (pi / 2)^(3 / 4)
normalization(sg) * 8 / √(π / 2)

f_1 = (T::Float64) -> GTO.F(1, T)
f_1 = function (T::Float64)
    T^2.0
end
function cubing(T::Float64)
    (T -> T^3.0)(T)
end
using BenchmarkTools
@benchmark f_1(1.0)
@benchmark cubing(1.0)
@code_lowered f_1(1.0)
@code_lowered cubing(1.0)

@benchmark f_1(1.0)
@benchmark 1.0^3
@time GTO.F(1, 1)

f = eval(:((α::Float64, a::Float64, b::Float64, c::Float64, T::Float64) -> @r_unroll(1, 1, 1, 0, α, a, b, c, T)))

function f_unrolled(α::Float64, a::Float64, b::Float64, c::Float64, T::Float64)
    @r_unroll(1, 1, 1, 0, α, a, b, c, T)
end
@time f_unrolled(1.0, 1.0, 1.0, 1.0, √(3))
