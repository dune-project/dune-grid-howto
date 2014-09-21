// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "unitcube.hh"
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

#include <dune/grid/onedgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {
    // make a grid
    UnitCube<Dune::OneDGrid,1> uc0;

    UnitCube<Dune::YaspGrid<3>,1> uc1;
    UnitCube<Dune::YaspGrid<2>,1> uc2;
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
#if ALBERTA_DIM==3
    UnitCube<Dune::AlbertaGrid<3,3>,1> uc9;
#endif
#endif
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
    UnitCube<Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>,1> uc8;
    UnitCube<Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>,1> uc10;
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
