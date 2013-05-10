// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/io/file/vtk/vtkwriter.hh> // VTK output routines
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

#include "unitcube.hh"
#include "functors.hh"
#include "integrateentity.hh"


//! adaptive refinement test
template<class Grid, class Functor>
void adaptiveintegration (Grid& grid, const Functor& f)
{
  // get grid view type for leaf grid part
  typedef typename Grid::LeafGridView GridView;
  // get iterator type
  typedef typename GridView::template Codim<0>::Iterator ElementLeafIterator;

  // get grid view on leaf part
  GridView gridView = grid.leafView();

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
    for (ElementLeafIterator it = gridView.template begin<0>();
         it!=gridView.template end<0>(); ++it)
      value += integrateEntity(*it,f,highorder);           /*@\label{aic:int1}@*/

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
      double lowresult=integrateEntity(*it,f,loworder);
      double highresult=integrateEntity(*it,f,highorder);
      double error = std::abs(lowresult-highresult);

      // max over whole grid
      maxerror = std::max(maxerror,error);

      // error on father entity
      double fatherlowresult=integrateEntity(*(it->father()),f,loworder);
      double fatherhighresult=integrateEntity(*(it->father()),f,highorder);
      double fathererror = std::abs(fatherlowresult-fatherhighresult);

      // local extrapolation
      double extrapolatederror = error*error/(fathererror+1E-30);
      maxextrapolatederror = std::max(maxextrapolatederror,extrapolatederror);
    }
    double kappa = std::min(maxextrapolatederror,0.5*maxerror);       /*@\label{aic:kappa1}@*/

    // mark elements for refinement
    for (ElementLeafIterator it = gridView.template begin<0>();       /*@\label{aic:mark0}@*/
         it!=gridView.template end<0>(); ++it)
    {
      double lowresult=integrateEntity(*it,f,loworder);
      double highresult=integrateEntity(*it,f,highorder);
      double error = std::abs(lowresult-highresult);
      if (error>kappa) grid.mark(1,*it);
    }                                                  /*@\label{aic:mark1}@*/

    // adapt the mesh
    grid.preAdapt();                                   /*@\label{aic:ref0}@*/
    grid.adapt();
    grid.postAdapt();                                  /*@\label{aic:ref1}@*/
  }

  // write grid in VTK format
  Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(gridView);
  vtkwriter.write( "adaptivegrid", Dune::VTK::appendedraw );
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

    // the GridSelector :: GridType is defined in gridtype.hh and is
    // set during compilation
    typedef GridSelector :: GridType Grid;

    // use unitcube from grids
    std::stringstream dgfFileName;
    dgfFileName << DUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH
      << "unitcube" << Grid::dimension << ".dgf";

    // create grid pointer
    GridPtr<Grid> gridPtr( dgfFileName.str() );

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
