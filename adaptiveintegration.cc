// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/io/file/vtk/vtkwriter.hh> // VTK output routines
#include <dune/common/mpihelper.hh> // include mpi helper class

// checks for defined gridtype and inlcudes appropriate dgfparser implementation
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>


#include "unitcube.hh"
#include "functors.hh"
#include "integrateentity.hh"


//! adaptive refinement test
template<class Grid, class Functor>
void adaptiveintegration (Grid& grid, const Functor& f)
{
  // get iterator type
  typedef typename Grid::template Codim<0>::LeafIterator ElementLeafIterator;

  // algorithm parameters
  const double tol=1E-8;
  const int loworder=1;
  const int highorder=3;

  // loop over grid sequence
  double oldvalue=1E100;
  for (int k=0; k<100; k++)
  {
    // compute integral on current mesh
    double value=0;                                      /*@\label{aic:int0}@*/
    for (ElementLeafIterator it = grid.template leafbegin<0>();
         it!=grid.template leafend<0>(); ++it)
      value += integrateentity(it,f,highorder);           /*@\label{aic:int1}@*/

    // print result
    double estimated_error = std::abs(value-oldvalue);
    oldvalue=value;       // save value for next estimate
    std::cout << "elements="
              << std::setw(8) << std::right
              << grid.size(0)
              << " integral="
              << std::scientific << std::setprecision(8)
              << value
              << " error=" << estimated_error
              << std::endl;

    // check convergence
    if (estimated_error <= tol*value)                  /*@\label{aic:finish}@*/
      break;

    // refine grid globally in first step to ensure
    // that every element has a father
    if (k==0)
    {
      grid.globalRefine(1);                            /*@\label{aic:gr}@*/
      continue;
    }

    // compute threshold for subsequent refinement
    double maxerror=-1E100;                            /*@\label{aic:kappa0}@*/
    double maxextrapolatederror=-1E100;
    for (ElementLeafIterator it = grid.template leafbegin<0>();
         it!=grid.template leafend<0>(); ++it)
    {
      // error on this entity
      double lowresult=integrateentity(it,f,loworder);
      double highresult=integrateentity(it,f,highorder);
      double error = std::abs(lowresult-highresult);

      // max over whole grid
      maxerror = std::max(maxerror,error);

      // error on father entity
      double fatherlowresult=integrateentity(it->father(),f,loworder);
      double fatherhighresult=integrateentity(it->father(),f,highorder);
      double fathererror = std::abs(fatherlowresult-fatherhighresult);

      // local extrapolation
      double extrapolatederror = error*error/(fathererror+1E-30);
      maxextrapolatederror = std::max(maxextrapolatederror,extrapolatederror);
    }
    double kappa = std::min(maxextrapolatederror,0.5*maxerror);       /*@\label{aic:kappa1}@*/

    // mark elements for refinement
    for (ElementLeafIterator it = grid.template leafbegin<0>();       /*@\label{aic:mark0}@*/
         it!=grid.template leafend<0>(); ++it)
    {
      double lowresult=integrateentity(it,f,loworder);
      double highresult=integrateentity(it,f,highorder);
      double error = std::abs(lowresult-highresult);
      if (error>kappa) grid.mark(1,it);
    }                                                  /*@\label{aic:mark1}@*/

    // adapt the mesh
    grid.preAdapt();                                   /*@\label{aic:ref0}@*/
    grid.adapt();
    grid.postAdapt();                                  /*@\label{aic:ref1}@*/
  }

  // write grid in VTK format
  Dune::VTKWriter<Grid> vtkwriter(grid);
  vtkwriter.write("adaptivegrid",Dune::VTKOptions::binaryappended);
}

//! supply functor
template<class Grid>
void dowork (Grid& grid)
{
  adaptiveintegration(grid,Needle<typename Grid::ctype,Grid::dimension>());
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    using namespace Dune;

    // use unitcube from grids
    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    GridPtr<GridType> gridPtr( dgfFileName.str() );

    // do the adaptive integration
    // NOTE: for structured grids global refinement will be used
    dowork( *gridPtr );
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
