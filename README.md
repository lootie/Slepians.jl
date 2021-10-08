# Slepians.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lootie.github.io/Slepians.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://lootie.github.io/Slepians.jl/dev)

Slepians.jl is a package that solves the concentration problem of Slepian,
Landau and Pollak, numerically. Loosely speaking, we find a function that is
both finite in extent, and has a Fourier transform which lies in a certain
domain. In a single dimension, with discrete sampling, the discrete prolate 
spheroidal sequences solve this optimization problem. 


# Installation

Slepians.jl is unregistered and relies on unregistered packages.  To avoid
difficulties in which Julia does not know the relevant URLS, I have created a
registry bbkt-reg.jl which will tell your installation where to find this
package and its dependencies. Begin with adding the registry using 

```
pkg> registry add https://github.com/lootie/bbkt-reg.jl
```

then one can simply add the Slepians.jl package as

```
pkg> add Slepians
```

# Examples

The jupyter notebooks in the `Examples/` directory illustrate the functionality
of this package.

# Funding

This material is based upon work supported by the U.S. Department of Energy, 
Office of Science, Office of Basic Energy Sciences.
