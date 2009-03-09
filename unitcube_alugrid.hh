// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALU3DGRID_HH
#define UNITCUBE_ALU3DGRID_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

// ALU3dGrid and ALU2dGrid simplex specialization.
// Note: element type determined by type
template<>
class UnitCube<Dune::ALUSimplexGrid<3,3>,1>
  : public BasicUnitCube< 3 >
{
public:
  typedef Dune::ALUSimplexGrid<3,3> GridType;

private:
  GridType * grid_;

public:
  UnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< 3 >::insertVertices( factory );
    BasicUnitCube< 3 >::insertSimplices( factory );
    grid_ = factory.createGrid( );
  }

  ~UnitCube()
  {
    delete grid_;
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

// ALU2SimplexGrid 2d specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALUSimplexGrid<2,2>,1>
{
public:
  typedef Dune::ALUSimplexGrid<2,2> GridType;

  UnitCube () : filename("grids/2dsimplex.alu"), grid_(filename.c_str())
  {}

  GridType& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  GridType grid_;
};

// ALU3dGrid hexahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALUCubeGrid<3,3>,1>
  : public BasicUnitCube< 3 >
{
public:
  typedef Dune::ALUCubeGrid<3,3> GridType;

private:
  GridType * grid_;

public:
  UnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< 3 >::insertVertices( factory );
    BasicUnitCube< 3 >::insertCubes( factory );
    grid_ = factory.createGrid( );
  }

  ~UnitCube()
  {
    delete grid_;
  }

  GridType &grid ()
  {
    return *grid_;
  }
};
#endif

#endif
