// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>

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
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  /*
     UnitCube<Dune::OneDGrid<1,1>,1> uc0;
     UnitCube<Dune::YaspGrid<3,3>,1> uc1;
     UnitCube<Dune::YaspGrid<2,2>,1> uc2;
     UnitCube<Dune::SGrid<1,1>,1> uc3;
     UnitCube<Dune::SGrid<2,2>,1> uc4;
     UnitCube<Dune::SGrid<3,3>,1> uc5;
     #if HAVE_UG
     UnitCube<Dune::UGGrid<3,3>,2> uc6;
     #endif
     #if HAVE_ALBERTA
     #if DUNE_PROBLEM_DIM==2
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

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
