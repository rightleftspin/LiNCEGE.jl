# LINCEGE
Put Badges here

LINCEGE.jl is a Julia implementation of Numerical Linked Cluster Expansions (NLCEs) for generic finite or infinite lattices. Utilizing LINCEGE.jl, you can easily generate clusters and reduced lattice constants for any lattice in order to perform NLCE.

## Examples
### Square Lattice Site Expansion

### Square Lattice Square Expansion

## Try out LINCEGE.jl
Follow the documented tutorials for LINCEGE.jl below.

[`Stable`](https://LINCEGE.github.io/LINCEGE.jl/stable/): Documentation for the latest version of the code

## Contact and Citation
To contact the authors of this package, please email either Pranav Seetharaman at [pjseetha@uwaterloo.ca](mailto:pjseetha@uwaterloo.ca) or Professor Ehsan Khatami at [ehsan.khatami@sjsu.edu](mailto:ehsan.khatami@sjsu.edu).

If LINCEGE.jl helped you in your research, please cite us using the following citation:
```bibtex
```

Funded by the U.S. Department of Energy: Grant Number DE-SC0022311


This section only needs to be run the first time you clone the repository

julia --project=.

using Pkg

Pkg.instantiate()

exit()

From here on out, you can just run

julia --project=. examples/ex-1.jl

or any other file in the repository.
