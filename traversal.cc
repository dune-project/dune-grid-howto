// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

// C/C++ includes
#include <iostream>           // for standard I/O

// Dune includes
#include "config.h"           // file constructed by ./configure script
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class


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

  // type of the GridView used for traversal
  // every grid exports a LeafGridView and a LevelGridView
  typedef typename G :: LeafGridView LeafGridView;     /*@\label{tc:lfgv}@*/

  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();         /*@\label{tc:lfv}@*/

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LeafGridView::template Codim<0>::Iterator ElementLeafIterator; /*@\label{tc:ittype}@*/

  // iterate through all entities of codim 0 at the leaves
  int count = 0;
  for (ElementLeafIterator it = leafView.template begin<0>();      /*@\label{tc:forel}@*/
       it!=leafView.template end<0>(); ++it)
  {                                                    /*@\label{tc:forel0}@*/
    Dune::GeometryType gt = it->type();       /*@\label{tc:reftype}@*/
    std::cout << "visiting leaf " << gt                /*@\label{tc:print}@*/
              << " with first vertex at " << it->geometry().corner(0)
              << std::endl;
    count++;                                           /*@\label{tc:count}@*/
  }                                                    /*@\label{tc:forel1}@*/

  std::cout << "there are/is " << count << " leaf element(s)" << std::endl;

  // Leafwise traversal of codim dim
  std::cout << std::endl;
  std::cout << "*** Traverse codim " << dim << " leaves" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LeafGridView :: template Codim<dim>
  :: Iterator VertexLeafIterator;                  /*@\label{tc:vertit}@*/

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (VertexLeafIterator it = leafView.template begin<dim>(); /*@\label{tc:forve}@*/
       it!=leafView.template end<dim>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    std::cout << "visiting " << gt
              << " at " << it->geometry().corner(0)
              << std::endl;
    count++;
  }
  std::cout << "there are/is " << count << " leaf vertices(s)"
            << std::endl;

  // Levelwise traversal of codim 0
  std::cout << std::endl;
  std::cout << "*** Traverse codim 0 level-wise" << std::endl;

  // type of the GridView used for traversal
  // every grid exports a LeafGridView and a LevelGridView
  typedef typename G :: LevelGridView LevelGridView;   /*@\label{tc:level0}@*/

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LevelGridView :: template Codim<0>
  :: Iterator ElementLevelIterator;

  // iterate through all entities of codim 0 on the given level
  for (int level=0; level<=grid.maxLevel(); level++)
  {
    // get the instance of the LeafGridView
    LevelGridView levelView = grid.levelGridView(level);

    count = 0;
    for (ElementLevelIterator it = levelView.template begin<0>();
         it!=levelView.template end<0>(); ++it)
    {
      Dune::GeometryType gt = it->type();
      std::cout << "visiting " << gt
                << " with first vertex at " << it->geometry().corner(0)
                << std::endl;
      count++;
    }
    std::cout << "there are/is " << count << " element(s) on level "
              << level << std::endl;
    std::cout << std::endl;
  }                                                    /*@\label{tc:level1}@*/
}                                                      /*@\label{tc:tra1}@*/


int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

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
