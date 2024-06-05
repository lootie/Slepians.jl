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
date: Mar 23, 2023
bibliography: paper.bib

---

# Summary

In a series of six landmark papers of Slepian, Landau, and Pollak, a set of functions
are described which solve what is known as the ``concentration problem". The
concentration problem refers to the problem of functions having limited support in,
say, time, while simultaneously having a Fourier transform for which the bulk of
its nonzero mass in a small constrained interval. While in some special discrete
cases, it is possible to solve a simple, possibly tridiagonal or Toeplitz eigenvalue problem, in others it is necessary
to use numerical integration. This package produces function values where the problem
dimension is up to three. We rely on the work of Simons, Wang,
[@SimonsWang2011] and others, for the two dimensional case but to our knowledge the
three dimensional case is novel to this package. In addition, we implement the
missing data Slepian sequences of Chave [@chave] as well as the generalized Slepian
sequences of Bronez [@bronez1988] which solve the 2D concentration problem on an
unequally spaced spatial grid.

# Mathematical Preliminaries

The concentration problem of Slepian, Landau and Pollak, seeks to find bandlimited functions
$g(t)$, i.e. functions satisfying 

$$ g(t) = \frac{1}{2 \pi} \int_{-W}^{W} G(\omega) e^{i\omega t} d\omega $$

where $G(\omega)$ denotes the Fourier transform of the function $g(t)$, i.e.

$$ G(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i\omega t} dt $$

that maximize the ratio 

$$ \lambda = \frac{\int_{-T}^{T} g^2(t)dt}{\int_{-\infty}^{\infty} g^2(t)dt}. $$

$\lambda$ is referred to as the concentration. It was shown that these functions
satisfy the Fredholm integral equation of the second kind

$$ \int_{-T}^{T} \frac{\sin [W (t - s)]}{\pi (t-s)} g(s)ds = \lambda g(t) $$

for real $t$. 

## Discrete prolate spheroidal sequences

The discrete prolate spheroidal sequences [@Slepian1978], $v_n^{(k)}(N,W)$, are the
discrete analog to the above functions $g(t)$. As such, they are the eigenvector
solutions to the eigenvalue problem

$$\sum_{m = 0}^{N-1}\frac{\sin 2\pi W(n-m)}{\pi(n-m)}v_m^{(k)}(N,W) = \lambda_k(N,W)\cdot v_n^{(k)}(N,W)$$

where $N$ is the length of the sequence, $n = 0, \ldots, N-1$, $W$ is the bandwidth,
$K$ is the number of tapers, and $\lambda_k(N,W)$ is the eigenvalue corresponding to
the $k$th taper. 

The dpss's solve the problem of concentrating, in frequency, the most amount of mass
under the interval $(-W,W)$, while having only finite extent in time. The first $2NW$
tapers have eigenvalues close to one, whereby their magnitude drops off rapidly
thereafter. 

Of note is that the dpss's are never computed using the above equation because of
numerical instability. Grunbaum [@grunbaum1981] showed that the associated differential operator
commutes with a symmetric tridiagonal matrix $T$, which is

$$ T_{nn} = [(N-1-2n)/2]^2 \cos (2\pi W), \quad n = 0, \ldots, N-1$$

on the diagonal and 

$$ T_{n, n + 1} = T_{n, n-1} = (n + 1)(N - n - 1)/2, \quad n = 0, \ldots, N-2$$

on the sub and super diagonals. The matrix $T$ has the same eigenvectors as the
problem above and the eigenvalues can be found by substitution into the original
equation. The tridiagonal formulation reduces the computation to one which is soluble
in O($N$) time, and the solution tends to be more accurate.

Show a figure as an example.

Generate a Durrani-Chapman filter

## Generalized Prolate Spheroidal Sequences 
 
The missing-data multitaper method requires its own set of optimal tapers
\cite{T82,bronez88,chave2019multitaper}, the \emph{missing-data prolate spheroidal
sequences} or missing-data Slepian sequences. These sequences extend the notion of
maximally bandlimited orthogonal sequences, \cite{S78}, to those sampled on a grid
where there is missing-data.  Define the bandwidth $W$ for the desired spectral
window, and let the sequences be given on the length $N$ index set
$\{t_n\}_{n=0}^{N-1}$ with unit sampling except where there are gaps. 

The missing-data Slepian sequences solve the eigenvalue problem, equations (21)-(23) in [@Chave2019],
$$ \label{eq:mdslep} 
  \lambda_k v^{k}_{n} = \sum_{m=0}^{N-1} \frac{\sin 2\pi W (t_n - t_m) }{\pi (t_n - t_m)} v^{k}_{m}.$$
The $k$th missing-data Slepian sequence is denoted $v^{k}_{n}$ at time $t_n$, and
and $\lambda_k$ is its associated eigenvalue, sorted in decreasing order
$1>\lambda_0>\lambda_1>\ldots>\lambda_{N-1}>0$. Note that both of these explicitly
depend on $N$ and $W$, but for simplicity this is suppressed in the notation. When
there are no missing values, the sinc-like matrix in \eqref{eq:mdslep} is
Toeplitz and has a special form which allows for numerically effective eigenvalue
routines.  

By applying the nonuniform discrete Fourier transform to the generalized data taper
in \eqref{eq:mdslep}, one obtains the missing-data prolate spheroidal functions, or
\emph{missing-data Slepian functions}

$$ V^{k}(f) = \sum_{n = 0}^{N-1} v^{k}_{n} e^{-2 \pi i t_n f}. $$

The missing-data Slepian sequences form an orthonormal set on the grid
$\{t_n\}_{n=0}^{N-1}$; and the missing-data Slepian functions are orthonormal on
$(-1/2,1/2)$ and orthogonal on $(-W,W)$. 

Show a figure as an example




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
