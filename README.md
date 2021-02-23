[![Build Status](https://github.com/temken/libphysica/workflows/Build%20Status/badge.svg)](https://github.com/temken/libphysica/actions)
[![codecov](https://codecov.io/gh/temken/libphysica/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/libphysica)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# libphysica
Static C++ library collecting functions, variables, and classes for application in scientific codes.

- *Linear_Algebra*: Basic linear algebra functionality (vectors and matrices).
- *Natural_Units*: A simple implementation of natural units.
- *Numerics*: Various special functions and numerical algorithms (integration, interpolation, root_finding, multidimensional minimization).
- *Statistics*: Includes PDFs and CDFs for a number of distributions, sampling techniques, rudimentary data analysis, and Kernel density estimation.
- *Utilities*: Some useful functions, that don't fit anywhere else. (Progress bars, import and export of data from files,...)

## DEPENDENCIES

- The unit tests are set up with the [googletest](https://github.com/google/googletest) framework, which is downloaded and installed automatically via [CMake](https://cmake.org/).

## VERSION HISTORY

- 23.02.2021: Release of version 0.1.0

## AUTHORS & CONTACT

The author of this library is Timon Emken.

For questions, bug reports, or other suggestions, please contact [emken@chalmers.se](mailto:emken@chalmers.se).

## LICENSE

This project is licensed under the MIT License - see the LICENSE file.

## ACKNOWLEDGEMENTS

A number of functions in the *Numerics* module use implementations from this book:

- [*Numerical Recipes 3rd Edition -  The Art of Scientific Computing*](https://en.wikipedia.org/wiki/Numerical_Recipes)  
W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery   
Cambridge University Press, (2007)
