using Printf

function Base.show(io::IO, g::CartesianGaussian)
    exp_string = @sprintf("%.3f", -α(g))
    center_string = "[" * join(map(n->@sprintf("%.3f", n), c(g)), ",") * "]"
    print(io, "x^$(i(g)) y^$(j(g)) z^$(k(g)) exp($(exp_string) (r - $(center_string)))")
end

function Base.show(io::IO, cg::ContractedGaussian)
    for (c, g) in cg
        print(io, @sprintf("%.3f ", c))
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