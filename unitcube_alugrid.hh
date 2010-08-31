// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_ALUGRID_HH
#define UNITCUBE_ALUGRID_HH

#include "unitcube.hh"

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

// ALU3dGrid and ALU2dGrid simplex specialization.
// Note: element type determined by type
template<int dim>
class UnitCube<Dune::ALUSimplexGrid<dim,dim>,1>
{
public:
  typedef Dune::ALUSimplexGrid<dim,dim> GridType;

private:
  Dune::shared_ptr<GridType> grid_;

public:
  UnitCube ()
  {
    Dune::FieldVector<typename GridType::ctype,dim> lowerLeft(0);
    Dune::FieldVector<typename GridType::ctype,dim> upperRight(1);
    Dune::array<unsigned int,dim> elements;
    std::fill(elements.begin(), elements.end(), 1);

    grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

// ALU3dGrid hexahedra specialization. Note: element type determined by type
template<>
class UnitCube<Dune::ALUCubeGrid<3,3>,1>
{
public:
  typedef Dune::ALUCubeGrid<3,3> GridType;

private:
  Dune::shared_ptr<GridType> grid_;

public:
  UnitCube ()
  {
    Dune::FieldVector<GridType::ctype,3> lowerLeft(0);
    Dune::FieldVector<GridType::ctype,3> upperRight(1);
    Dune::array<unsigned int,3> elements = {1,1,1};

    grid_ = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
  }

  GridType &grid ()
  {
    return *grid_;
  }
};
#endif

#endif
