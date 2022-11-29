# GTO

## Introduction

This package aims to easy-to-use Gaussian integrals.  The code prioritizes
flexibility over performance since the targeted audience is method developers
including myself, for whom performance is usually a secondary concern compared
to flexibility.

The API of this package provides what a package user naively would want as
oppose to what gives the optimal performance. This does not mean that the
performance is completely sacrificed; it just means that it should be done
at the cost of complicating the API.

The theory implemented is [McMurchie & Davidson 1978](https://doi.org/10.1016/0021-9991\(78\)90092-X)
and [Boys 1949](https://royalsocietypublishing.org/doi/10.1098/rspa.1950.0036) 

## Things that need work

1. More extensive testing.
2. Implement a cache for better performance.
3. Full documentation (the readme seems sufficient).

## Usage

### Load a basis set

The basis set is downloaded from [BSE](https://www.basissetexchange.org/).
Only `.json` files are supported at the moment. 

To load a basis set, pass in the filename, which can be found in `basis_set_bundle-json-bib`.

```julia
b = load_basis("sto-3g.1.json")
```

### Atoms and AOs

To create the atom, pass in the atomic number and the coordinates.  The unit is
currently in Bohr by default, but I don't know if that is standard.
It seems that `pyscf` uses Angstrom, which can be converted through `BOHR_TO_ANGSTROM`.

```julia
a_1 = Atom(6, [0, 0, 0])
a_2 = Atom(6, [2 / BOHR_TO_ANGSTROM, 0, 0])
```

To create AOs for the atom, load them from the basis set `b`.
The result is `shells` of contracted Gaussians. The basis 
can then be constructed from the shells.

```julia
shells = make_gaussians(b, a_1, a_2)
basis = vcat(get_basis.(shells)...)
```

### Gaussians

#### Cartesian Gaussian

One typically does not need to construct this by hand, but it is good to know
what they are. A Cartesian Gaussian with angular momentum `i,j,k`, exponent `α`,
and center `c` can be constructed as

```julia
c = CartesianGaussian(0, 0, 0, 1, [1, 0, 0])
#                     i  j  k  α,  c
```

#### Hermite Gaussian

The Hermite Gaussians can be constructed similarly to the Cartesian Gaussians.
They are not used in GTO theories to represent the basis, but to represent products
of Cartesian Gaussians.

```julia
c = HermiteGaussian(0, 0, 0, 1, [1, 0, 0])
#                   i  j  k  α,  c
```

#### Spherical Gaussians

Not yet implemented. Don't need it for the moment.

#### Contracted Gaussians

These are just sums of Gaussians. Sums of Cartesian Gaussians are used to
represent contracted basis, whereas contracted Hermite Gaussians are used 
to represent products of (possibly contracted) Cartesian Gaussians.

TODO: Come up with an API.

#### The Product of two Gaussians

To multiply two (contracted) Cartesian Gaussians `p` and `q` together, use

```julia
gaussian_product(p, q)
p * q   # Syntax sugar.
```

The result will be a contracted Hermite Gaussian.
There is current no support (or need) for multiplying Hermite Gaussians.

### Integrals

The integrals necessary for HF and possibly mildly beyond are implemented.

#### Overlap

The overlap integrals are performed as

```julia
overlap_integral(p, q)
p' * q    # Syntax sugar.
S = [p' * q for p in basis, q in basis] # The overlap matrix.
```

#### Kinetic

The kinetic integrals are performed as 
```julia
kinetic_integral(p, q)
p' * ∇² * q     # Syntax sugar.
```

#### Nuclear

The nuclear integral for a few `Atom`s can be written as
```julia
nuclear_potential(p * q, a_1) + nuclear_potential(p * q, a_2)
p' * VNuc(a_1, a_2) * q     # Syntax sugar.
```

#### Two electrons 

Two electron integrals are performed as
```julia
two_electron_integral(p*q, r*s)
(p*q | r*s)   # Syntax sugar.
```