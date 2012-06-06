// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <dune/common/mpihelper.hh> // include mpi helper class
#include <dune/grid/io/file/dgfparser/dgfparser.hh>


#include "elementdata.hh"
#include "vertexdata.hh"
#include "functors.hh"
#include "unitcube.hh"


#ifdef GRIDDIM
const int dimGrid = GRIDDIM;
#else
const int dimGrid = 2;
#endif


//! supply functor
template<class Grid>
void dowork ( Grid &grid, int refSteps = 5 )
{
  // make function object
  Exp<typename Grid::ctype,Grid::dimension> f;

  // refine the grid
  grid.globalRefine( refSteps );

  // call the visualization functions
  elementdata(grid,f);
  vertexdata(grid,f);
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try
  {
    if( argc > 1 )
    {
      typedef Dune::GridSelector::GridType DGFGridType;
      // create grid pointer
      Dune :: GridPtr< DGFGridType > gridPtr( argv[ 1 ] );
      dowork( *gridPtr, 3 );
    }

    /*
       UnitCube<Dune::OneDGrid,1> uc0;
       UnitCube<Dune::YaspGrid<dimGrid>,1> uc1;

       #if HAVE_UG
       UnitCube< Dune::UGGrid< dimGrid >, 2 > uc2;
       dowork( uc2.grid(), 3 );
       #endif

       #if HAVE_ALBERTA
       {
       UnitCube< Dune::AlbertaGrid< dimGrid, dimGrid >, 1 > unitcube;
       // note: The 3d cube cannot be bisected recursively
       dowork( unitcube.grid(), (dimGrid < 3 ? 6 : 0) );
       }
       #endif // #if HAVE_ALBERTA
     */

    UnitCube< Dune::SGrid< dimGrid, dimGrid >, 1 > uc4;
    dowork( uc4.grid(), 3 );

#if HAVE_ALUGRID
    UnitCube< Dune::ALUGrid< dimGrid, dimGrid, Dune::simplex,
            Dune::nonconforming > , 1 > uc5;
    dowork( uc5.grid(), 3 );

#if GRIDDIM == 2 || GRIDDIM == 3
    UnitCube< Dune::ALUGrid< dimGrid, dimGrid, Dune::cube,
            Dune::nonconforming > , 1 > uc6;
    dowork( uc6.grid(), 3 );
#endif // #if GRIDDIM == 2 || GRIDDIM == 3
#endif // #if HAVE_ALUGRID
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
