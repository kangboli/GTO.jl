export Cache, cache, cache!


# abstract type AbstractAtom end


mutable struct Cache
    n_atoms::Int
    n_cartesian::Int
    cartesian_normalization_factors::Vector{Float64}
    gaussian_atom::Vector{<:AbstractAtom}
end

Cache() = Cache(0, 0, Vector{Float64}(), Vector{AbstractAtom}())
generate_atom_cache_number!(c::Cache) = c.n_atoms += 1
genreate_cartesian_cache_number!(c::Cache) = c.n_cartesian += 1

# The grand global cache.
ggc = Cache()

# mutable struct Atom <: AbstractAtom
#     atomic_number::Int
#     coordinates::Vector{Float64}
#     cache_number::Int
# end

# function Atom(atomic_number::Int, coordinates::Vector{<:Real})
#     Atom(atomic_number, float.(coordinates), -1)
# end

function cache!(cache::Cache, a::Atom)
    a.cache_number = generate_atom_cache_number!(cache)
    return a
end

