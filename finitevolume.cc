// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

#include "vtkout.hh"
#include "transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"

//===============================================================
// the time loop function working for all types of grids
//===============================================================

template<class G>
void timeloop (const G& grid, double tend)
{
  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,Dune::MCMGElementLayout>
  mapper(grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid,mapper,c);                           /*@\label{fvc:init}@*/
  vtkout(grid,c,"concentration",0,0.0);

  // now do the time steps
  double t=0,dt;
  int k=0;
  const double saveInterval = 0.1;
  double saveStep = 0.1;
  int counter = 1;

  while (t<tend)                                       /*@\label{fvc:loop0}@*/
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
      vtkout(grid,c,"concentration",counter,t);     /*@\label{fvc:file}@*/

      // increase counter and saveStep for next interval
      saveStep += saveInterval;
      ++counter;
    }

    // print info about time, timestep size and counter
    std::cout << "s=" << grid.size(0)
              << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
  }                                                    /*@\label{fvc:loop1}@*/
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

    // use unitcube from dgf grids
    std::stringstream dgfFileName;
    dgfFileName << DUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH
      << "unitcube" << Grid::dimension << ".dgf";

    // create grid pointer
    GridPtr<Grid> gridPtr( dgfFileName.str() );

    // grid reference
    Grid& grid = *gridPtr;

    // half grid width 4 times
    int level = 4 * DGFGridInfo<Grid>::refineStepsForHalf();

    // refine grid until upper limit of level
    grid.globalRefine(level);

    // do time loop until end time 0.5
    timeloop(grid, 0.5);
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
