// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include "vtkout.hh"

#include "unitcube.hh"
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
  bool contains (int codim, Dune::GeometryType gt)
  {
    if (codim==0) return true;
    return false;
  }
};

template<class G>
void timeloop (const G& grid, double tend)
{
  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
  mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);
  vtkout(grid,c,"concentration",0);

  // now do the time steps
  double t=0,dt;
  int k=0;
  const int modulo=5;
  while (t<tend)
  {
    k++;
    evolve(grid,mapper,c,t,dt);
    if (k%modulo==0) vtkout(grid,c,"concentration",k/modulo);
    std::cout << "k=" << k << " t=" << t << " dt=" << dt << std::endl;
    t += dt;
  }

  // output results
  vtkout(grid,c,"concentration",k/modulo);
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv)
{
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  UnitCube<Dune::YaspGrid<2,2>,1> uc;
  UnitCube<Dune::OneDGrid<1,1>,1> uc0;
  UnitCube<Dune::SGrid<1,1>,1> uc1;
#if HAVE_UG
  UnitCube<Dune::UGGrid<2,2>,2> uc2;
#endif
#if HAVE_ALBERTA
#if DUNE_PROBLEM_DIM==2
  UnitCube<Dune::AlbertaGrid<2,2>,1> uc3;
#endif
#endif

  uc0.grid().globalRefine(7);
  timeloop(uc0.grid(),0.5);

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
