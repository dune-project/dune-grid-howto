// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// Dune includes
#include "config.h"           // file constructed by ./configure script
#include "dune/grid/sgrid.hh" // load sgrid definition
#include "dune/grid/common/gridinfo.hh" // definition of gridinfo

int main()
{
  // make a grid
  const int dim=3;
  typedef Dune::SGrid<dim,dim> GridType;
  Dune::FieldVector<int,dim> N(3);
  Dune::FieldVector<GridType::ctype,dim> L(-1.0);
  Dune::FieldVector<GridType::ctype,dim> H(1.0);
  GridType grid(N,L,H);

  // print some information about the grid
  Dune::gridinfo(grid);

  // done
  return 0;
}
