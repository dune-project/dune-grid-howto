// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/common/mpihelper.hh> // include mpi helper class

#include "vtkout.hh"
#include "unitcube.hh"
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"
#include "finitevolumeadapt.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

template<class G>
void timeloop (G& grid, double tend, int lmin, int lmax)
{
  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,Dune::MCMGElementLayout>
  mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);
  for (int i=grid.maxLevel(); i<lmax; i++)
  {
    if (grid.maxLevel()>=lmax) break;
    finitevolumeadapt(grid,mapper,c,lmin,lmax,0);           /*@\label{afv:in}@*/
    initialize(grid,mapper,c);
  }

  // write initial data
  vtkout(grid,c,"concentration",0,0);

  // variables for time, timestep etc.
  double dt, t=0;
  double saveStep = 0.1;
  const double saveInterval = 0.1;
  int counter = 0;
  int k = 0;

  std::cout << "s=" << grid.size(0) << " k=" << k << " t=" << t << std::endl;
  while (t<tend)
  {
    // augment time step counter
    ++k;

    // apply finite volume scheme
    evolve(grid,mapper,c,t,dt);

    // augment time
    t += dt;

    // check if data should be written
    if (t >= saveStep)
    {
      // write data
      vtkout(grid,c,"concentration",counter,t);

      // increase counter and saveStep for next interval
      saveStep += saveInterval;
      ++counter;
    }

    // print info about time, timestep size and counter
    std::cout << "s=" << grid.size(0)
              << " k=" << k << " t=" << t << " dt=" << dt << std::endl;

    // for unstructured grids call adaptation algorithm
    finitevolumeadapt(grid,mapper,c,lmin,lmax,k);        /*@\label{afv:ad}@*/
  }

  // write last time step
  vtkout(grid,c,"concentration",counter,tend);

  // write
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv)
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

    // grid reference
    Grid& grid = *gridPtr;

    // minimal allowed level during refinement
    int minLevel = 2 * DGFGridInfo<Grid>::refineStepsForHalf();

    // refine grid until upper limit of level
    grid.globalRefine(minLevel);

    // maximal allowed level during refinement
    int maxLevel = minLevel + 3 * DGFGridInfo<Grid>::refineStepsForHalf();

    // do time loop until end time 0.5
    timeloop(grid, 0.5, minLevel, maxLevel);
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
