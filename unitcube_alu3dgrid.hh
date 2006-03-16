// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALU3DGRID_HH
#define UNITCUBE_ALU3DGRID_HH

#if HAVE_ALUGRID
#include <dune/grid/alu3dgrid.hh>

// ALU3dGrid tetrahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALU3dGrid<3,3,Dune::tetra>,1>
{
public:
  typedef Dune::ALU3dGrid<3,3,Dune::tetra> GridType;

  UnitCube () : filename("grids/cube.tetra"), grid_(filename.c_str())
  {}

  Dune::ALU3dGrid<3,3,Dune::tetra>& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  Dune::ALU3dGrid<3,3,Dune::tetra> grid_;
};

// ALU3dGrid hexahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALU3dGrid<3,3,Dune::hexa>,1>
{
public:
  typedef Dune::ALU3dGrid<3,3,Dune::hexa> GridType;

  UnitCube () : filename("grids/cube.hexa"), grid_(filename.c_str())
  {}

  Dune::ALU3dGrid<3,3,Dune::hexa>& grid ()
  {
    return grid_;
  }

private:
  std::string filename;
  Dune::ALU3dGrid<3,3,Dune::hexa> grid_;
};
#endif

#endif
