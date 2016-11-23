DUNE-grid-howto
===============

DUNE, the Distributed and Unified Numerics Environment is a modular toolbox
for solving partial differential equations with grid-based methods.

The main intention is to create slim interfaces allowing an efficient use of
legacy and/or new libraries. Using C++ techniques DUNE allows to use very
different implementation of the same concept (i.e. grid, solver, ...) under
a common interface with a very low overhead.

DUNE was designed with flexibility in mind. It supports easy discretization
using methods, like Finite Elements, Finite Volume and also Finite
Differences. Through separation of data structures DUNE allows fast Linear
Algebra like provided in the ISTL module, or usage of external libraries
like blas.

This package contains the basic DUNE how-to classes and the grid-howto
document. It is an introduction to the usage of the DUNE grid interface.
Some techniques might be a little bit outdated.

More information
----------------

Check dune-common for more details concerning dependencies, known bugs,
license and installation.
