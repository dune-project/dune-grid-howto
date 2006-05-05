// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include "dune/grid/io/file/vtk/vtkwriter.hh"
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
    double value=0;
    for (ElementLeafIterator it = grid.template leafbegin<0>();
         it!=grid.template leafend<0>(); ++it)
      value += integrateentity(it,f,highorder);

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
    if (estimated_error <= tol*value)
      break;

    // refine grid globally in first step to ensure
    // that every element has a father
    if (k==0)
    {
      grid.globalRefine(1);
      continue;
    }

    // compute threshold for subsequent refinement
    double maxerror=-1E100;
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
      double fathererror = std::abs(lowresult-highresult);

      // local extrapolation
      double extrapolatederror = error*error/(fathererror+1E-30);
      maxextrapolatederror = std::max(maxextrapolatederror,extrapolatederror);
    }
    double kappa = std::min(maxextrapolatederror,0.5*maxerror);

    // mark elements for refinement
    for (ElementLeafIterator it = grid.template leafbegin<0>();
         it!=grid.template leafend<0>(); ++it)
    {
      double lowresult=integrateentity(it,f,loworder);
      double highresult=integrateentity(it,f,highorder);
      double error = std::abs(lowresult-highresult);
      if (error>kappa) grid.mark(1,it);
    }

    // adapt the mesh
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
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
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

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

#ifdef HAVE_UG
  dowork(uc6.grid());
#else
  dowork(uc2.grid());
#endif

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
