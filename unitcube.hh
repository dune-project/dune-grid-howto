// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UNITCUBE_HH
#define UNITCUBE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// default implementation for any template parameter
template<typename T, int variant>                      /*@\label{uc:uc0}@*/
class UnitCube
{
public:
  typedef T GridType;

  // constructor throwing exception
  UnitCube ()
  {
    DUNE_THROW(Dune::Exception,"no specialization for this grid available");
  }

  T& grid ()
  {
    return grid_;
  }

private:
  // the constructed grid object
  T grid_;
};                                                     /*@\label{uc:uc1}@*/

// include basic unitcube using GridFactory concept
#include "basicunitcube.hh"

// include specializations
#include "unitcube_onedgrid.hh"
#include "unitcube_sgrid.hh"
#include "unitcube_yaspgrid.hh"
#include "unitcube_uggrid.hh"
#include "unitcube_albertagrid.hh"
#include "unitcube_alugrid.hh"

#endif
