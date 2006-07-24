// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALU3DGRID_HH
#define UNITCUBE_ALU3DGRID_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>

// ALU3dGrid tetrahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALUSimplexGrid<3,3>,1>
{
public:
  typedef Dune::ALUSimplexGrid<3,3> GridType;

  UnitCube () : filename("grids/cube.tetra"), grid_(filename.c_str())
  {}

  GridType& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  GridType grid_;
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
{
public:
  typedef Dune::ALUCubeGrid<3,3> GridType;

  UnitCube () : filename("grids/cube.hexa"), grid_(filename.c_str())
  {}

  GridType& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  GridType grid_;
};
#endif

#endif
