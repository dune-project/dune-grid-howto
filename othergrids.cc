// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include "config.h"
#include "dune/grid/common/gridinfo.hh"
#include "unitcube.hh"

template<typename G>
void test (const G& grid)
{
  typedef typename G::ctype ct;

  // assignability of leaf iterator
  typedef typename G::template Codim<0>::LeafIterator LeafIterator;
  LeafIterator leafi(grid.template leafbegin<0>());
  LeafIterator leafj(grid.template leafend<0>());
  leafj = leafi;

  // equality comparison
  if (leafi==leafj)
    std::cout << "equal" << std::endl;

  // assignability of level iterator
  typedef typename G::template Codim<0>::LevelIterator LevelIterator;
  LevelIterator leveli(grid.template lbegin<0>(0));
  LevelIterator levelj(grid.template lend<0>(0));
  levelj = leveli;

  // equality comparison
  if (leveli==levelj)
    std::cout << "equal" << std::endl;

  // assignability of intersection iterator
  typedef typename G::template Codim<0>::IntersectionIterator IntersectionIterator;
  IntersectionIterator interi(leafi->ibegin());
  IntersectionIterator interj(leafi->iend());
  interj = interi;

  // equality comparison
  if (interi==interj)
    std::cout << "equal" << std::endl;

  // assignability of hierarchic iterator
  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
  HierarchicIterator hi(leafi->hbegin(100));
  HierarchicIterator hj(leafi->hend(100));
  hj = hi;

  // equality comparison
  if (hi==hj)
    std::cout << "equal" << std::endl;

  // entity pointer
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
  EntityPointer epi(leafi->father());
  EntityPointer epj(leveli);
  EntityPointer epk(leafi);
  EntityPointer epm(hi);
  epj = epi;
  epk = leafi;

  // equality comparison
  if (epi==epj)
    std::cout << "equal" << std::endl;

  //   if (leafi->boundaryId());
}

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
  UnitCube<Dune::UGGrid<3,3>,2> uc6;
#endif
#if HAVE_ALBERTA
  UnitCube<Dune::AlbertaGrid<3,3>,1> uc7;
#endif
#if HAVE_ALUGRID
  UnitCube<Dune::ALU3dGrid<3,3,Dune::hexa>,1> uc8;
#endif

  // print some information about the grid
  test(uc1.grid());

#if HAVE_MPI
  MPI_Finalize();
#endif

  // done
  return 0;
}
