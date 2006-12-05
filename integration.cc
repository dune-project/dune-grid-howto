// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// Dune includes
#include "config.h"
#include <dune/grid/sgrid.hh>

#include "unitcube.hh"
#include "functors.hh"
#include "integrateentity.hh"

//! uniform refinement test
template<class Grid>
void uniformintegration (Grid& grid)
{
  // function to integrate
  Exp<typename Grid::ctype,Grid::dimension> f;

  // get iterator type
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;

  // loop over grid sequence
  double oldvalue=1E100;
  for (int k=0; k<10; k++)
  {
    // compute integral with some order
    double value = 0.0;
    LeafIterator eendit = grid.template leafend<0>();
    for (LeafIterator it = grid.template leafbegin<0>(); it!=eendit; ++it)
      value += integrateentity(it,f,1);                /*@\label{ic:call}@*/

    // print result and error estimate
    std::cout << "elements="
              << std::setw(8) << std::right
              << grid.size(0)
              << " integral="
              << std::scientific << std::setprecision(12)
              << value
              << " error=" << std::abs(value-oldvalue)
              << std::endl;

    // save value of integral
    oldvalue=value;

    // refine all elements
    grid.globalRefine(1);
  }
}

int main(int argc, char **argv)
{
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // start try/catch block to get error messages from dune
  try {
    // make a grid
    UnitCube<Dune::OneDGrid<1,1>,1> uc0;
    UnitCube<Dune::SGrid<2,2>,1> uc1;
    UnitCube<Dune::YaspGrid<2,2>,1> uc2;
    UnitCube<Dune::YaspGrid<3,3>,1> uc3;

    // integrate and compute error with extrapolation
    uniformintegration(uc2.grid());
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
