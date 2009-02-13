// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_UGGRID_HH
#define UNITCUBE_UGGRID_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>

// UGGrid 3d, variant 1 (hexahedra) specialization
template<int variant>
class UnitCube<Dune::UGGrid<3>,variant>
{
public:
  typedef Dune::UGGrid<3> GridType;

  UnitCube ()
  {
    //   Start grid creation
    Dune::GridFactory<GridType> factory;

    //   Insert vertices
    Dune::FieldVector<double,3> pos;

    pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);

    if (variant==1) {

      // Insert element
      std::vector<unsigned int> cornerIDs(8);
      for (int i=0; i<8; i++)
        cornerIDs[i] = i;

      factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3), cornerIDs);

    } else {

      // Insert elements
      std::vector<unsigned int> cornerIDs(4);

      cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

      cornerIDs[0] = 1;  cornerIDs[1] = 3;  cornerIDs[2] = 2;  cornerIDs[3] = 7;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

      cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 4;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

      cornerIDs[0] = 1;  cornerIDs[1] = 7;  cornerIDs[2] = 4;  cornerIDs[3] = 5;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

      cornerIDs[0] = 4;  cornerIDs[1] = 7;  cornerIDs[2] = 2;  cornerIDs[3] = 6;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

    }

    //   Finish initialization
    grid_ = factory.createGrid();
  }

  GridType& grid ()
  {
    return *grid_;
  }

private:
  GridType* grid_;
};

// UGGrid 2d,
template<int variant>
class UnitCube<Dune::UGGrid<2>, variant>
{
public:
  typedef Dune::UGGrid<2> GridType;

  UnitCube ()
  {
    //   Start grid creation
    Dune::GridFactory<GridType> factory;

    //   Insert vertices
    Dune::FieldVector<double,2> pos;

    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);

    pos[0] = 1;  pos[1] = 0;
    factory.insertVertex(pos);

    pos[0] = 0;  pos[1] = 1;
    factory.insertVertex(pos);

    pos[0] = 1;  pos[1] = 1;
    factory.insertVertex(pos);

    if (variant==1) {

      // Insert element
      std::vector<unsigned int> cornerIDs(4);
      cornerIDs[0] = 0;
      cornerIDs[1] = 1;
      cornerIDs[2] = 2;
      cornerIDs[3] = 3;

      factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);

    } else {

      // Insert element
      std::vector<unsigned int> cornerIDs(3);

      cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);

      cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);

    }

    //   Finish initialization
    grid_ = factory.createGrid();
  }

  GridType& grid ()
  {
    return *grid_;
  }

private:
  GridType* grid_;
};

#endif

#endif
