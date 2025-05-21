using Dates
using JSON

export load_basis,
    make_gaussians,
    get_basis,
    BasisIO,
    generate_basis


struct ElectronShell
    function_type::String
    angular_momentum::Vector{Int}
    region::String
    exponents::Vector{Float64}
    coefficients::Matrix{Float64}
end

struct Reference
    reference_keys::Vector{String}
    reference_description::String
end

struct BasisElement
    electron_shells::Vector{ElectronShell}
    references::Vector{Reference}
end

struct BasisIO
    revision_data::Date
    family::String
    role::String
    revision_description::String
    elements::Dict{Int,BasisElement}
    names::Vector{String}
    name::String
    auxiliaries::Dict{String,<:Any}
    molssi_bse_schema::Dict{String,<:Any}
    version::String
    description::String
    function_types::Vector{String}
    tags::Vector{<:Any}
end


function load_basis(filename::String)
    basis_set_dict = JSON.parsefile(joinpath(pathof(GTO)[1:end-11],
        "basis_set_bundle-json-bib", filename))

    elements = Dict(
        parse(Int, k) => BasisElement(
            [
                ElectronShell(
                    s["function_type"],
                    s["angular_momentum"],
                    s["region"],
                    map(f -> parse(Float64, f), s["exponents"]),
                    (f -> parse(Float64, f)).(hcat(s["coefficients"]...)),
                ) for s in v["electron_shells"]
            ],
            [
                Reference(r["reference_keys"], r["reference_description"]) for
                r in v["references"]
            ],
        ) for (k, v) in basis_set_dict["elements"]
    )

    BasisIO(
        Date(basis_set_dict["revision_date"]),
        basis_set_dict["family"],
        basis_set_dict["role"],
        basis_set_dict["revision_description"],
        elements,
        basis_set_dict["names"],
        basis_set_dict["name"],
        basis_set_dict["auxiliaries"],
        basis_set_dict["molssi_bse_schema"],
        basis_set_dict["version"],
        basis_set_dict["description"],
        basis_set_dict["function_types"],
        basis_set_dict["tags"],
    )
end

function make_gaussians(basis::BasisIO, atoms::Vararg{<:AbstractAtom})

    function make_gaussians_single(atom)
        element = basis.elements[atomic_number(atom)]

        function make_shell(shell::ElectronShell)
            if shell.function_type == "gto" || shell.function_type == "gto_cartesian"
                Dict(l => make_angular(shell, l, atom) for l in shell.angular_momentum)
            elseif shell.function_type == "gto_spherical"
                Dict(l => make_spherical(shell, l, atom) for l in shell.angular_momentum)
            else
                error("function type $(shell.function_type) is not supported.")
            end
        end

        make_shell.(element.electron_shells)
    end

    vcat(map(make_gaussians_single, atoms)...)
end

make_shells = make_gaussians

function make_angular(shell::ElectronShell, l::Int, atom::AbstractAtom)
    result = []
    for (i, j, k) in collect(Iterators.product(0:l, 0:l, 0:l))
        i + j + k == l || continue
        cartesians = map(α -> CartesianGaussian(i, j, k, α, coordinates(atom)), shell.exponents)
        # We make the choice to carried the normalization in the contraction coefficients.
        g = contract(shell.coefficients[:, findfirst(i -> i == l, shell.angular_momentum)] .* normalization.(cartesians), cartesians, true)
        push!(result, g)
    end
    return result
end

function make_spherical(shell::ElectronShell, l::Int, atom::AbstractAtom)
    result = []
    #= for m in -l:l
        sphericals = map(α -> SphericalGaussian(l, m, l, α, coordinates(atom)), shell.exponents)
        g = contract(shell.coefficients[:, findfirst(i -> i == l, shell.angular_momentum)] .* normalization.(sphericals), sphericals)
        push!(result, g)
    end =#

    spherical_sym = map(α -> SphericalGaussian(l, 0, l, α, coordinates(atom)), shell.exponents)
    coeffs = shell.coefficients[:, findfirst(i -> i == l, shell.angular_momentum)] .* normalization.(spherical_sym)
    push!(result, contract(coeffs, spherical_sym))

    #  1/√2 (Y(l,m) ± Y(l,-m)) to make the basis real
    for m in 1:l
        postives = map(α -> SphericalGaussian(l, m, l, α, coordinates(atom)), shell.exponents)
        netagives = map(α -> SphericalGaussian(l, -m, l, α, coordinates(atom)), shell.exponents)
        #= coeffs = shell.coefficients[:, findfirst(i -> i == l, shell.angular_momentum)] .* normalization.(postives) =#
        g_plus = contract(repeat(coeffs, 2) ./ sqrt(2), vcat(postives, netagives))
        g_minus = contract(vcat(coeffs, -coeffs) ./ sqrt(-(Complex(2))), vcat(postives, netagives))
        push!(result, g_plus)
        push!(result, g_minus)
    end
    return result
end

function get_basis(shell::Dict{Int,<:Any})
    vcat(collect(values(shell))...)
end

function generate_basis(basis_set::BasisIO, molecule::Vararg{AbstractAtom})
    basis = []
    for a in molecule
        shells = make_shells(basis_set, a)

        gaussian_layer(l) = vcat(map(sh -> get(sh, l, []), shells)...)
        gaussians = vcat(map(gaussian_layer, 0:7)...)
        append!(basis, gaussians)
    end
    return basis
end
