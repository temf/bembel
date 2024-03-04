<!-- This file is part of Bembel, the higher order C++ boundary element library.

Copyright (C) 2024 see <http://www.bembel.eu>

It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
source code is subject to the GNU General Public License version 3 and
provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
information. -->
# Bembel
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/temf/bembel)
[![DOI](https://zenodo.org/badge/173278911.svg)](https://zenodo.org/badge/latestdoi/173278911)
[![Documentation](https://img.shields.io/badge/Documentation-Doxygen-blue)](https://temf.github.io/bembel/Doxy_out/html/index.html)
[![coverage report](https://temf.github.io/bembel/coverage/report/coverage.svg)](https://temf.github.io/bembel/coverage/report/report.html)
## Table of contents
1. [Introduction](#introduction)
2. [What is a Bembel?](#whatis)
3. [Features](#features)
4. [Documentation](#doc)
5. [Known Bugs and Upcoming Features](#bugs)
6. [Contributing to Bembel](#contributing)
7. [Publications, Preprints, and how to cite](#publications)

## 1. Introduction <a name="introduction"></a>

Bembel is the *Boundary Element Method Based Engineering Library* written in C++ library featuring higher order isogeometric Galerkin boundary element methods for Laplace, Helmholtz, and Maxwell problems.
Bembel is compatible with geometries from the Octave NURBS package, and provides an interface to the Eigen template library for linear algebra operations.
For computational efficiency, it applies an embedded fast multipole method tailored to the isogeometric analysis framework and a parallel matrix assembly based on OpenMP.

## 2. What is a Bembel?<a name="whatis"></a>

A traditional Hessian ceramic, as depicted in our logo. 
Quoting [Wikipedia](https://en.wikipedia.org/wiki/Apfelwein):

> *Most establishments also serve Apfelwein by the Bembel (a specific Apfelwein
> jug), much like how beer can be purchased by the pitcher in many countries. The
> paunchy Bembel (made from salt-glazed stoneware) usually has a basic grey colour
> with blue-painted detailing.*

## 3. Features <a name="features"></a>

Current key features include

* Header-only implementation,
* Already implemented Laplace, Helmholtz and Maxwell single layer operator,
easily expendable to other operators,
* Arbitrary parametric mappings for the geometry representation, by default
realized as NURBS-mappings from files generated by the
[NURBS package](https://octave.sourceforge.io/nurbs/),
* Arbitrary-order B-Spline functions as Ansatz spaces, as in
the framework of isogeometric analysis for electromagnetics,
* An embedded interpolation-based fast multipole method for compression,
equivalent to the H2 matrix format,
* openMP parallelized matrix assembly,
* Full compatibility with the [Eigen](http://eigen.tuxfamily.org/) linear algebra library.

## 4. Documentation <a name="doc"></a>

A good place to start are the examples in the `examples/` folder, together with publication [[1]](#1). Apart from that, a [Doxygen documentation](https://temf.github.io/bembel/Doxy_out/html/index.html) is available.

## 5. Known Bugs and Upcoming Features <a name="bugs"></a>

For a list of known bugs and upcoming features, please have a look at 
the issue tracker on github.

## 6. Contributing to Bembel <a name="contributing"></a>

Any contribution to this project in fixing a bug or implementing a feature is welcome.
Create a fork of this repository and create a pull request after your finished the implementation of the feature.
To successfully merge your pull request you should follow our [Coding Guidelines](https://temf.github.io/bembel/Doxy_out/html/codingGuidelines.html)


## 7. How to cite <a name="publications"></a>

<a name="1">[1]</a> J. Dölz, H. Harbrecht, S. Kurz, M. Multerer, S. Schöps, and F. Wolf. *Bembel: The Fast Isogeometric Boundary Element C++ Library for Laplace, Helmholtz, and Electric Wave Equation*. In: SoftwareX, 11, 10476.
[![doi badge](https://img.shields.io/badge/DOI-10.1016/j.softx.2020.100476-blue)](https://doi.org/10.1016/j.softx.2020.100476)

```bib
@article{Bembel2020,
doi = {https://doi.org/10.1016/j.softx.2020.100476},
year = {2020}, volume = {11}, pages = {100476},
issn = {2352-7110},
author = {J. Dölz and H. Harbrecht and S. Kurz and M. Multerer and S. Schöps and F. Wolf},
title = {{Bembel}: The fast isogeometric boundary element {C++} library for {Laplace}, {Helmholtz}, and electric wave equation},
journal = {Software {X}}}
```