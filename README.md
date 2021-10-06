# Slepians.jl



# Installation

Slepians.jl relies on KernelMatrices.jl which is unregistered, however, one can use the custom registry bbkt-reg.jl which will tell your installation where to find KernelMatrices.jl. 

```
pkg> registry add https://github.com/lootie/bbkt-reg.jl
```

then one can simply add the Slepians.jl package as

```
pkg> add Slepians
```

or, if you prefer

```
pkg> add https://github.com/lootie/Slepians.jl
```

