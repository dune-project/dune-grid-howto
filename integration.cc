// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// Dune includes
#include "config.h"           // file constructed by ./configure script
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

#include "functors.hh"
#include "integrateentity.hh"

//! uniform refinement test
template<class Grid>
void uniformintegration (Grid& grid)
{
  // function to integrate
  Exp<typename Grid::ctype,Grid::dimension> f;

  // get GridView on leaf grid - type
  typedef typename Grid :: LeafGridView GridView;

  // get GridView instance
  GridView gridView = grid.leafView();

  // get iterator type
  typedef typename GridView :: template Codim<0> :: Iterator LeafIterator;

  // loop over grid sequence
  double oldvalue=1E100;
  for (int k=0; k<10; k++)
  {
    // compute integral with some order
    double value = 0.0;
    LeafIterator eendit = gridView.template end<0>();
    for (LeafIterator it = gridView.template begin<0>(); it!=eendit; ++it)
      value += integrateEntity(*it,f,1);                /*@\label{ic:call}@*/

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
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    using namespace Dune;

    // the GridSelector :: GridType is defined in gridtype.hh and is
    // set during compilation
    typedef GridSelector :: GridType Grid;

    // use unitcube from grids
    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube" << Grid :: dimension << ".dgf";

    // create grid pointer
    GridPtr<Grid> gridPtr( dgfFileName.str() );

    // integrate and compute error with extrapolation
    uniformintegration( *gridPtr );
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

  // done
  return 0;
}
