// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALBERTAGRID_HH
#define UNITCUBE_ALBERTAGRID_HH

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>

// AlbertaGrid 2d, variant 1 (2 triangles) specialization
#if ALBERTA_DIM == 2 && ALBERTA_WORLD_DIM == 2
template<>
class UnitCube<Dune::AlbertaGrid<2,2>,1>
{
public:
  typedef Dune::AlbertaGrid<2,2> GridType;

  UnitCube () : grid_("grids/2dgrid.al")
  {}

  Dune::AlbertaGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::AlbertaGrid<2,2> grid_;
};
#endif

// AlbertaGrid 3d, variant 1 (6 tetrahedra) specialization
#if ALBERTA_DIM == 3 && ALBERTA_WORLD_DIM == 3
template<>
class UnitCube<Dune::AlbertaGrid<3,3>,1>
{
public:
  typedef Dune::AlbertaGrid<3,3> GridType;

  UnitCube () : grid_("grids/3dgrid.al")
  {}

  Dune::AlbertaGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::AlbertaGrid<3,3> grid_;
};
#endif
#endif
#endif
