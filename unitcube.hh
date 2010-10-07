// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_HH
#define UNITCUBE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

// default implementation for any template parameter
template<typename T, int variant>                      /*@\label{uc:uc0}@*/
class UnitCube
{
public:
  typedef T GridType;

  static const int dim = GridType::dimension;

  // constructor throwing exception
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
                  << variant << " of unit cube not implemented." );
    }
  }

  T& grid ()
  {
    return *grid_;
  }

private:
  // the constructed grid object
  Dune::shared_ptr<T> grid_;
};                                                     /*@\label{uc:uc1}@*/


// include specializations
#include "unitcube_sgrid.hh"
#include "unitcube_yaspgrid.hh"
#include "unitcube_albertagrid.hh"
#include "unitcube_alugrid.hh"

#endif
