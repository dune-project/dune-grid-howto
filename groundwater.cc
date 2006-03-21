// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include "unitcube.hh"
#include "groundwaterproblem.hh"
#include "groundwateradapt.hh"

int main(int argc, char **argv)
{
#if HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // make a grid
  UnitCube<Dune::OneDGrid<1,1>,1> uc0;
  UnitCube<Dune::YaspGrid<3,3>,1> uc1;
  UnitCube<Dune::YaspGrid<2,2>,1> uc2;
  UnitCube<Dune::SGrid<1,1>,1> uc3;
  UnitCube<Dune::SGrid<2,2>,1> uc4;
  UnitCube<Dune::SGrid<3,3>,1> uc5;
#if HAVE_UG
  UnitCube<Dune::UGGrid<2,2>,2> uc7b;
#endif
#if HAVE_ALBERTA
#if DUNE_PROBLEM_DIM==2
  UnitCube<Dune::AlbertaGrid<2,2>,1> uc8;
#endif
#if DUNE_PROBLEM_DIM==3
  UnitCube<Dune::AlbertaGrid<3,3>,1> uc9;
#endif
#endif
#if HAVE_ALUGRID
  UnitCube<Dune::ALU3dGrid<3,3,Dune::hexa>,1> uc10;
#endif

  // do the work
  groundwateradapt(uc7b.grid(),20,1,false,true);

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
