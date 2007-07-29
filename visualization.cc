// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <dune/common/mpihelper.hh> // include mpi helper class


#include "elementdata.hh"
#include "vertexdata.hh"
#include "functors.hh"
#include "unitcube.hh"

//! supply functor
template<class Grid>
void dowork (Grid& grid)
{
  // make function object
  Exp<typename Grid::ctype,Grid::dimension> f;

  // refine the grid
  grid.globalRefine(5);

  // call the visualization functions
  elementdata(grid,f);
  vertexdata(grid,f);
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    /*
       UnitCube<Dune::OneDGrid,1> uc0;
       UnitCube<Dune::YaspGrid<3,3>,1> uc1;
       UnitCube<Dune::YaspGrid<2,2>,1> uc2;
       UnitCube<Dune::SGrid<1,1>,1> uc3;
       UnitCube<Dune::SGrid<2,2>,1> uc4;
       UnitCube<Dune::SGrid<3,3>,1> uc5;
       #if HAVE_UG
       UnitCube<Dune::UGGrid<3>,2> uc6;
       #endif
       #if HAVE_ALBERTA
       #if ALBERTA_DIM==2
       UnitCube<Dune::AlbertaGrid<2,2>,1> uc7;
       #endif
       #endif
     */
    UnitCube<Dune::SGrid<2,2>,1> uc4;
    dowork(uc4.grid());

#if HAVE_ALUGRID
    UnitCube<Dune::ALUSimplexGrid<3,3> ,1> uc8;
    dowork(uc8.grid());
#endif
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
