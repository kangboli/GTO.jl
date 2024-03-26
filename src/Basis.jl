using Dates
using JSON

export load_basis,
make_gaussians,
get_basis,
BOHR_TO_ANGSTROM,
charge,
BasisIO

const BOHR_TO_ANGSTROM = 0.529177248994098

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
    vcat(map(a->make_gaussians(basis, a), atoms)...)
end

function make_gaussians(basis::BasisIO, atom::AbstractAtom)
    element = basis.elements[atomic_number(atom)]
    function make_shell(shell::ElectronShell)
        Dict(l => make_angular(shell, l) for l in shell.angular_momentum)
    end

    function make_angular(shell::ElectronShell, l::Int)
        map(filter(ls -> sum(ls) == l, collect(product(0:l, 0:l, 0:l)))) do (i, j, k)
            cartesians = map(α -> CartesianGaussian(i, j, k, α, coordinates(atom)), shell.exponents)
            # We make the choice to carried the normalization in the contraction coefficients.
            contract(shell.coefficients[:, l+1] .* normalization.(cartesians), cartesians , true)
        end
    end
    make_shell.(element.electron_shells)
end


function get_basis(shell::Dict{Int, <:Any})
    vcat(collect(values(shell))...)
end

