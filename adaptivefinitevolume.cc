// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include "dune/grid/common/mcmgmapper.hh" // mapper class

#include "unitcube.hh"
#include "transportproblem.hh"
#include "initialize.hh"
#include "evolve.hh"
#include "finitevolumeadapt.hh"
#include "vtkconcentration.hh"

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
void timeloop (G& grid, double tend, int lmin, int lmax)
{
  // get leaf index set type needed for mapper
  typedef typename G::template Codim<0>::LeafIndexSet IS;

  // make a mapper for codim 0 entities in the leaf grid
  Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,P0Layout>
  mapper(grid,grid.leafIndexSet());

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);
  for (int i=grid.maxLevel(); i<lmax; i++)
  {
    finitevolumeadapt(grid,mapper,c,lmin,lmax);
    initialize(grid,mapper,c);
  }

  double dt, t=0;
  int k=0;
  while (t<tend)
  {
    k++;
    evolve(grid,mapper,c,t,dt);
    t += dt;
    std::cout << "s=" << grid.size(0) << " k=" << k
              << " t=" << t << " dt=" << dt << std::endl;
    finitevolumeadapt(grid,mapper,c,lmin,lmax);
  }

  // output results
  vtkconcentration(grid,c,0);
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
#if HAVE_UG
  UnitCube<Dune::UGGrid<2,2>,2> uc2;
#endif
#if HAVE_ALBERTA
  UnitCube<Dune::AlbertaGrid<2,2>,1> uc3;
#endif
  uc3.grid().globalRefine(8);
  timeloop(uc3.grid(),0.5,8,17);
  //     uc2.grid().globalRefine(5);
  //     timeloop(uc2.grid(),0.5,5,9);

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
