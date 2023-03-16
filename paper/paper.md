---
title: 'Slepians.jl: A Julia package to generate optimal functions which are 
simultaneously constrained in space and spectral domains'
tags:
  - Julia
  - David Slepian
  - bandlimited
  - timelimited
  - concentration problem
  - eigenvalue problem
  - spectral analysis
  - reproducing kernel
  - Spherical harmonics
  - Sturm-Liouville problem
authors:
  - name: Charlotte L. Haley
    orcid: 0000-0003-3996-773X
    affiliation: "1"
    email: "haley@anl.gov" 
affiliations:
 - name: Argonne National Laboratory
   index: 1
date: Mar 16, 2023
bibliography: paper.bib

---

# Summary

In a series of six landmark papers of Slepian, Landau, and Pollak, a set of functions
are described which solve what is known as the ``concentration problem". The
concentration problem refers to the problem of functions having limited support in,
say, time, while simultaneously having the bulk of its nonzero mass in a small
constrained interval centered at the origin. While in some special discrete cases, it
is possible to solve a simple eigenvalue problem, in others it is necessary to use
numerical integration. This package produces function values where the problem is
from one to three dimensional. We rely on the work of Simons, Wang, [@SimonsWang2011]
and others, for the two dimensional case but to our knowledge the three dimensional
case is novel to this package. In addition, we implement the missing data Slepian
sequences of Chave [@chave] as well as the generalized Slepian sequences of Bronez
[@bronez1988] which solve the 2D concentration problem on an unequally spaced spatial 
grid.

# Statement of need

`Slepians.jl` is a Julia package for producing simultaneously time and bandlimited
sequences, as well as space and spectrally limited functions up to three dimensions.
These functions were introduced by Slepian, Landau, and Pollak in the 1960s and have
had particular application in diverse fields, perhaps most notably as foundational
for the multitaper method for spectrum analysis [@T82] which extends both to 2D
Cartesian domains [@SimonsWang2011] as well as the sphere [@simons2006]. 

The high-level character of Julia allows for widely readable and extendible codes,
while the low-level functionality provides speed and efficiency. The `Multitaper.jl`
package provides a user-friendly implementation of many of the basic concepts.
Implementations of higher-dimensional Slepian tapers on Cartesian domains
[@SimonsWang2011] [@Geoga2018]. In addition, we provide tutorial-style notebooks to
allow accessibility to those new to these concepts or to Julia in general.

`Slepians.jl` has been used in the context of 3D pair distribution function analysis
of crystal structure. 

# Other software

The Matlab software of Frederik Simons's research group is available on github under
\url{https://github.com/csdms-contrib}, of particular note is [@slepian_foxtrot]
which relies on [@slepian_alpha] and [@slepian_delta]. These pacakages

# To contribute

We welcome input of any kind via bitbucket issues or by pull requests.
Support requests can kindly be directed to haley@anl.gov.

# Acknowledgements

We acknowledge contributions from Mihai Anitescu, David J. Thomson, Sally
Dodson-Robinson during the writing of these codes. We also gratefully acknowledge the
help of our reviewers in editing the code repository.

This work was supported by the U.S. Department of Energy, Office of Science, Advanced
Scientific Computing Research, under contract number DE-AC02-06CH11357.

The submitted manuscript has been created by UChicago Argonne, LLC, Operator of
Argonne National Laboratory (“Argonne”). Argonne, a U.S. Department of Energy Office
of Science laboratory, is operated under Contract No. DE-AC02-06CH11357. The U.S.
Government retains for itself, and others acting on its behalf, a paid-up
nonexclusive, irrevocable worldwide license in said article to reproduce, prepare
derivative works, distribute copies to the public, and perform publicly and display
publicly, by or on behalf of the Government. The Department of Energy will provide
public access to these results of federally sponsored research in accordance with the
DOE Public Access Plan. http://energy.gov/downloads/doe-public-access-plan

# References
