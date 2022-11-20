export Atom, Cache, atomic_number, cache, coordinates, cache_number


abstract type AbstractAtom end

atomic_number(a::AbstractAtom) = a.atomic_number
cache(a::AbstractAtom) = a.cache
cache_number(a::AbstractAtom) = a.cache_number
coordinates(a::AbstractAtom) = a.coordinates

mutable struct Cache
    n_atoms::Int
end

Cache() = Cache(0)
generate_cache_number!(c::Cache) = c.n_atoms += 1


struct Atom <: AbstractAtom
    atomic_number::Int
    coordinates::Vector{Float64}
    cache::Union{Cache}
    cache_number::Int
end

function Atom(atomic_number::Int, coordinates::Vector{<:Real}, cache=Cache())
    Atom(atomic_number, float.(coordinates), cache, generate_cache_number!(cache))
end

