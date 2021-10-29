[![Build Status](https://github.com/temken/libphysica/workflows/Build%20Status/badge.svg)](https://github.com/temken/libphysica/actions)
[![codecov](https://codecov.io/gh/temken/libphysica/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/libphysica)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# libphysica
[![DOI](https://zenodo.org/badge/261012574.svg)](https://zenodo.org/badge/latestdoi/261012574)

Static C++ library collecting functions, variables, and classes for application in scientific codes.

- *Integration*: Numerical integration / quadrature (for both one- and multi-dimensional integrals).
- *Linear_Algebra*: Basic linear algebra functionality (vectors and matrices).
- *Natural_Units*: A simple implementation of natural units.
- *Numerics*: Various special functions and numerical algorithms (interpolation, root finding, multidimensional minimization).
- *Special_Functions*: Implementation of various functions.
- *Statistics*: Includes PDFs and CDFs for a number of distributions, sampling techniques, rudimentary data analysis, and Kernel density estimation.
- *Utilities*: Some useful functions, that don't fit anywhere else. (Progress bars, import and export of data from files,...)

## DEPENDENCIES

- The unit tests are set up with the [googletest](https://github.com/google/googletest) framework, which is downloaded and installed automatically via [CMake](https://cmake.org/).
- The numerical integration methods rely partially on [boost](https://www.boost.org/).
- The *Configuration* class uses the library [libconfig](https://hyperrealm.github.io/libconfig/).

## INSTALLATION

An example on how to include *libphysica* into your CMake build can be found in [this repository](https://github.com/temken/template_cpp_cmake_libphysica), a C++ template code which is built with CMake and automatically downloads and includes this library during the build.

## VERSION HISTORY

- 23.02.2021: Release of version 0.1.0

## AUTHORS & CONTACT

The author of this library is Timon Emken.

For questions, bug reports, or other suggestions, please contact [timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se).

## LICENSE

This project is licensed under the MIT License - see the LICENSE file.

## ACKNOWLEDGEMENTS

A number of functions use implementations from this book:

- [*Numerical Recipes 3rd Edition -  The Art of Scientific Computing*](https://en.wikipedia.org/wiki/Numerical_Recipes)  
W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery   
Cambridge University Press, (2007)
