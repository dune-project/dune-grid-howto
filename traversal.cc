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
template<class G>                                      /*@\label{tc:tra0}@*/
void traversal (G& grid)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;                        /*@\label{tc:dim}@*/

  // type used for coordinates in the grid
  // such a type is exported by every grid implementation
  typedef typename G::ctype ct;                        /*@\label{tc:ct}@*/

  // Leaf Traversal
  std::cout << "*** Traverse codim 0 leaves" << std::endl;

  // the grid has an iterator providing the access to
  // all elements (better codim 0 entities) which are leafs
  // of the refinement tree.
  // Note the use of the typename keyword and the traits class
  typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator; /*@\label{tc:ittype}@*/

  // iterate through all entities of codim 0 at the leafs
  int count = 0;
  for (ElementLeafIterator it = grid.template leafbegin<0>(); /*@\label{tc:forel}@*/
       it!=grid.template leafend<0>(); ++it)
  {                                                    /*@\label{tc:forel0}@*/
    Dune::GeometryType gt = it->geometry().type();       /*@\label{tc:reftype}@*/
    std::cout << "visiting leaf " << gt                /*@\label{tc:print}@*/
              << " with first vertex at " << it->geometry()[0]
              << std::endl;
    count++;                                           /*@\label{tc:count}@*/
  }                                                    /*@\label{tc:forel1}@*/

  std::cout << "there are/is " << count << " leaf element(s)" << std::endl;

  // Leafwise traversal of codim dim
  std::cout << std::endl;
  std::cout << "*** Traverse codim " << dim << " leaves" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename G::template Codim<dim>::LeafIterator VertexLeafIterator; /*@\label{tc:vertit}@*/

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (VertexLeafIterator it = grid.template leafbegin<dim>(); /*@\label{tc:forve}@*/
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
  typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator; /*@\label{tc:level0}@*/

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
  }                                                    /*@\label{tc:level1}@*/
}                                                      /*@\label{tc:tra1}@*/


int main()
{
  // start try/catch block to get error messages from dune
  try {
    // make a grid
    const int dim=2;
    typedef Dune::SGrid<dim,dim> GridType;
    Dune::FieldVector<int,dim> N(1);
    Dune::FieldVector<GridType::ctype,dim> L(-1.0);
    Dune::FieldVector<GridType::ctype,dim> H(1.0);
    GridType grid(N,L,H);

    // refine all elements once using the standard refinement rule
    grid.globalRefine(1);                                /*@\label{tc:refine}@*/

    // traverse the grid and print some info
    traversal(grid);                                     /*@\label{tc:call}@*/
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
