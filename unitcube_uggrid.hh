// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_UGGRID_HH
#define UNITCUBE_UGGRID_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/amirameshreader.hh>

// UGGrid 3d, variant 1 (hexahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3,3>,1>
{
public:
  typedef Dune::UGGrid<3,3> GridType;

  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3,3> >::read(grid_,"grids/ug3dhexagrid.am");
  }

  Dune::UGGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3,3> grid_;
};

// UGGrid 3d, variant 2 (tetrahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3,3>,2>
{
public:
  typedef Dune::UGGrid<3,3> GridType;

  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3,3> >::read(grid_,"grids/ug3dtetragrid.am");
  }

  Dune::UGGrid<3,3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3,3> grid_;
};

// UGGrid 2d, variant 1 (quadrilaterals) specialization
template<>
class UnitCube<Dune::UGGrid<2,2>,1>
{
public:
  typedef Dune::UGGrid<2,2> GridType;

  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<2,2> >::read(grid_,"grids/quadgrid.am");
  }

  Dune::UGGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2,2> grid_;
};

// UGGrid 2d, variant 2 (triangles) specialization
template<>
class UnitCube<Dune::UGGrid<2,2>,2>
{
public:
  typedef Dune::UGGrid<2,2> GridType;

  UnitCube () : grid_(800,10)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<2,2> >::read(grid_,"grids/trianggrid.am");
  }

  Dune::UGGrid<2,2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2,2> grid_;
};
#endif

#endif
