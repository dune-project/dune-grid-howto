// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_UGGRID_HH
#define UNITCUBE_UGGRID_HH

#include "unitcube.hh"

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

template< int dim, int variant >
class UnitCube< Dune::UGGrid< dim >, variant >
{
public:
  typedef Dune::UGGrid< dim > GridType;

private:
  Dune::shared_ptr<GridType> grid_;

public:
  UnitCube ()
  {
    Dune::FieldVector<typename GridType::ctype,dim> lowerLeft(0);
    Dune::FieldVector<typename GridType::ctype,dim> upperRight(1);
    Dune::array<unsigned int,dim> elements;
    std::fill(elements.begin(), elements.end(), 1);

    switch (variant) {
    case 1 :
      grid_ = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft, upperRight, elements);
      break;
    case 2 :
      grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft, upperRight, elements);
      break;
    default :
      DUNE_THROW( Dune::NotImplemented, "Variant "
                  << variant << " of UG unit cube not implemented." );
    }
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

#endif // #if HAVE_UG

#endif
