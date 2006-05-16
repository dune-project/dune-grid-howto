// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// C/C++ includes
#include <iostream>           // for standard I/O

// Dune includes
#include "config.h"           // file constructed by ./configure script
#include <dune/grid/sgrid.hh> // load sgrid definition

// example for a generic algorithm that traverses
// the entities of a given mesh in various ways
template<class G>
void traversal (G& grid)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;

  // type used for coordinates in the grid
  // such a type is exported by every grid implementation
  typedef typename G::ctype ct;

  // Leaf Traversal
  std::cout << "*** Traverse codim 0 leaves" << std::endl;

  // the grid has an iterator providing the access to
  // all elements (better codim 0 entities) which are leafs
  // of the refinement tree.
  // Note the use of the typename keyword and the traits class
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;

  // iterate through all entities of codim 0 at the leafs
  int count = 0;
  for (ElementLeafIterator it = grid.template leafbegin<0>();
       it!=grid.template leafend<0>(); ++it)
  {
    Dune::GeometryType gt = it->geometry().type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << it->geometry()[0]
              << std::endl;
    count++;
  }

  std::cout << "there are/is " << count << " leaf element(s)" << std::endl;

  // Leafwise traversal of codim dim
  std::cout << std::endl;
  std::cout << "*** Traverse codim " << dim << " leaves" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename G::template Codim<dim>::LeafIterator VertexLeafIterator;

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (VertexLeafIterator it = grid.template leafbegin<dim>();
       it!=grid.template leafend<dim>(); ++it)
  {
    Dune::GeometryType gt = it->geometry().type();
    std::cout << "visiting " << gt
              << " at " << it->geometry()[0]
              << std::endl;
    count++;
  }
  std::cout << "there are/is " << count << " leaf vertices(s)"
            << std::endl;

  // Levelwise traversal of codim 0
  std::cout << std::endl;
  std::cout << "*** Traverse codim 0 level-wise" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator;

  // iterate through all entities of codim 0 on the given level
  for (int level=0; level<=grid.maxLevel(); level++)
  {
    count = 0;
    for (ElementLevelIterator it = grid.template lbegin<0>(level);
         it!=grid.template lend<0>(level); ++it)
    {
      Dune::GeometryType gt = it->geometry().type();
      std::cout << "visiting " << gt
                << " with first vertex at " << it->geometry()[0]
                << std::endl;
      count++;
    }
    std::cout << "there are/is " << count << " element(s) on level "
              << level << std::endl;
    std::cout << std::endl;
  }
}


int main()
{
  // make a grid
  const int dim=2;
  typedef Dune::SGrid<dim,dim> GridType;
  Dune::FieldVector<int,dim> N(1);
  Dune::FieldVector<GridType::ctype,dim> L(-1.0);
  Dune::FieldVector<GridType::ctype,dim> H(1.0);
  GridType grid(N,L,H);

  // refine all elements once using the standard refinement rule
  grid.globalRefine(1);

  // traverse the grid and print some info
  traversal(grid);

  // done
  return 0;
}
