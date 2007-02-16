// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_UGGRID_HH
#define UNITCUBE_UGGRID_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>

#if HAVE_AMIRAMESH
#include <dune/grid/io/file/amirameshreader.hh>
// UGGrid 3d, variant 1 (hexahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3>,1>
{
public:
  typedef Dune::UGGrid<3> GridType;

  UnitCube () : grid_(800)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3> >::read(grid_,"grids/ug3dhexagrid.am");
  }

  Dune::UGGrid<3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3> grid_;
};

// UGGrid 3d, variant 2 (tetrahedra) specialization
template<>
class UnitCube<Dune::UGGrid<3>,2>
{
public:
  typedef Dune::UGGrid<3> GridType;

  UnitCube () : grid_(800)
  {
    Dune::AmiraMeshReader<Dune::UGGrid<3> >::read(grid_,"grids/ug3dtetragrid.am");
  }

  Dune::UGGrid<3>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<3> grid_;
};
#endif

// UGGrid 2d, variant 1 (quadrilaterals) specialization
template<>
class UnitCube<Dune::UGGrid<2>,1>
{
public:
  typedef Dune::UGGrid<2> GridType;

  UnitCube () : grid_(800)
  {
    //   Start grid creation
    grid_.createBegin();

    //   Insert vertices
    Dune::FieldVector<double,2> pos;

    pos[0] = 0;  pos[1] = 0;
    grid_.insertVertex(pos);

    pos[0] = 1;  pos[1] = 0;
    grid_.insertVertex(pos);

    pos[0] = 0;  pos[1] = 1;
    grid_.insertVertex(pos);

    pos[0] = 1;  pos[1] = 1;
    grid_.insertVertex(pos);

    // Insert element
    std::vector<unsigned int> cornerIDs(4);
    cornerIDs[0] = 0;
    cornerIDs[1] = 1;
    cornerIDs[2] = 2;
    cornerIDs[3] = 3;

    grid_.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);

    //   Finish initialization
    grid_.createEnd();
  }

  Dune::UGGrid<2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2> grid_;
};

// UGGrid 2d, variant 2 (triangles) specialization
template<>
class UnitCube<Dune::UGGrid<2>,2>
{
public:
  typedef Dune::UGGrid<2> GridType;

  UnitCube () : grid_(800)
  {
    //   Start grid creation
    grid_.createBegin();

    //   Insert vertices
    Dune::FieldVector<double,2> pos;

    pos[0] = 0;  pos[1] = 0;
    grid_.insertVertex(pos);

    pos[0] = 1;  pos[1] = 0;
    grid_.insertVertex(pos);

    pos[0] = 0;  pos[1] = 1;
    grid_.insertVertex(pos);

    pos[0] = 1;  pos[1] = 1;
    grid_.insertVertex(pos);

    // Insert element
    std::vector<unsigned int> cornerIDs(3);

    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    grid_.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);

    cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
    grid_.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);

    //   Finish initialization
    grid_.createEnd();
  }

  Dune::UGGrid<2>& grid ()
  {
    return grid_;
  }

private:
  Dune::UGGrid<2> grid_;
};
#endif

#endif
