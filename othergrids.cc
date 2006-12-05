// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include "unitcube.hh"

int main(int argc, char **argv)
{
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // start try/catch block to get error messages from dune
  try {
    // make a grid
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
#if DUNE_PROBLEM_DIM==3
    UnitCube<Dune::AlbertaGrid<3,3>,1> uc9;
#endif
#endif
#if HAVE_ALUGRID
    UnitCube<Dune::ALUCubeGrid<3,3>,1> uc8;
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

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
