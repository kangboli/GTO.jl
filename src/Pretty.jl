using Printf

function Base.show(io::IO, g::CartesianGaussian)
    exp_string = @sprintf("%.3f", -α(g))
    center_string = "[" * join(map(n -> @sprintf("%.3f", n), c(g)), ",") * "]"
    print(io, "x^$(i(g)) y^$(j(g)) z^$(k(g)) exp($(exp_string) |r - $(center_string)|²)")
end

function Base.show(io::IO, g::SphericalGaussian)
    exp_string = @sprintf("%.3f", -α(g))
    center_string = "[" * join(map(n -> @sprintf("%.3f", n), c(g)), ",") * "]"
    print(io, "Y$(get_l(g))$(get_m(g))(r) exp($(exp_string) |r - $(center_string)|²)")
end

function Base.show(io::IO, cg::ContractedGaussian)
    for (c, g) in cg
        if isa(c, Float64)
            print(io, @sprintf("%.3f ", c))
        else
            print(io, @sprintf("(%.3f + %.3f im) ", real(c), imag(c)))
        end
        show(io, g)
        println(io)
    end
end

function Base.show(io::IO, ::MIME"text/plain", cgs::Vector{CCG})
    for cg in cgs
        show(io, cg)
        println(io)
    end
end
